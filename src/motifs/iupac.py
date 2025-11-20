"""
IUPAC code handling for DNA and protein motifs.

This module provides dictionaries mapping IUPAC codes to the set of concrete
letters they represent, and a helper that converts an IUPAC-encoded motif
string into a regular expression pattern.
"""

from __future__ import annotations

from typing import Dict

# Unambiguous + ambiguous DNA IUPAC codes.
# Values are *strings* listing the concrete bases represented by each code.
DNA_IUPAC: Dict[str, str] = {
    # Unambiguous bases
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "U": "U",
    # Ambiguous bases
    "R": "AG",      # purine (A or G)
    "Y": "CT",      # pyrimidine (C or T)
    "S": "GC",      # strong (G or C)
    "W": "AT",      # weak (A or T)
    "K": "GT",      # keto (G or T)
    "M": "AC",      # amino (A or C)
    "B": "CGT",     # not A
    "D": "AGT",     # not C
    "H": "ACT",     # not G
    "V": "ACG",     # not T
    "N": "ACGT",    # any base
}

# Standard 20 amino acids + common ambiguous protein codes.
PROTEIN_IUPAC: Dict[str, str] = {
    # Unambiguous amino acids
    "A": "A",
    "C": "C",
    "D": "D",
    "E": "E",
    "F": "F",
    "G": "G",
    "H": "H",
    "I": "I",
    "K": "K",
    "L": "L",
    "M": "M",
    "N": "N",
    "P": "P",
    "Q": "Q",
    "R": "R",
    "S": "S",
    "T": "T",
    "V": "V",
    "W": "W",
    "Y": "Y",
    # Ambiguous / special codes
    "B": "DN",      # D or N
    "Z": "EQ",      # E or Q
    "X": "ACDEFGHIKLMNPQRSTVWY",  # any amino acid
}


def _lookup_table(kind: str) -> Dict[str, str]:
    """Return the appropriate IUPAC table for ``kind``.

    Parameters
    ----------
    kind:
        Either ``"dna"`` or ``"protein"`` (case-insensitive).

    Raises
    ------
    ValueError
        If an unknown kind is supplied.
    """
    kind_lower = kind.lower()
    if kind_lower == "dna":
        return DNA_IUPAC
    if kind_lower == "protein":
        return PROTEIN_IUPAC
    raise ValueError(f"Unknown IUPAC kind: {kind!r} (expected 'dna' or 'protein')")


def iupac_to_regex(pattern: str, kind: str) -> str:
    """Translate an IUPAC motif string into a regular expression.

    This function is intentionally *conservative*: it only expands known
    IUPAC codes to character classes. Characters that are not recognised
    as IUPAC codes (including regex metacharacters like ``[`` or ``*``)
    are copied through untouched.

    Examples
    --------
    >>> iupac_to_regex("ATNG", "dna")
    'AT[ACGT]G'
    >>> iupac_to_regex("N[^P][ST][^P]", "protein")
    'N[^P][ST][^P]'
    """
    table = _lookup_table(kind)
    out_parts: list[str] = []

    META_CHARS = ".^$*+?{}[]|()\\\\"  # 기본 regex 메타 문자들

    for ch in pattern:
        # Leave regex metacharacters alone – we only expand bare IUPAC codes.
        if ch in META_CHARS:
            out_parts.append(ch)
            continue

        code = table.get(ch.upper())
        if code is None:
            # Not an IUPAC code for this alphabet, just keep it
            out_parts.append(ch)
        elif len(code) == 1:
            # Unambiguous code – use the literal letter
            out_parts.append(code)
        else:
            # Ambiguous code – expand to a character class
            out_parts.append("[" + code + "]")

    return "".join(out_parts)


__all__ = ["DNA_IUPAC", "PROTEIN_IUPAC", "iupac_to_regex"]
