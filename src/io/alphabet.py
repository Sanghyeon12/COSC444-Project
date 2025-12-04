from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Set

AlphabetType = Literal["DNA", "PROTEIN", "UNKNOWN"]


DNA_ALPHABET: Set[str] = {"A", "C", "G", "T"}

# 20 standard amino acids
PROTEIN_ALPHABET: Set[str] = {
    "A", "R", "N", "D", "C",
    "Q", "E", "G", "H", "I",
    "L", "K", "M", "F", "P",
    "S", "T", "W", "Y", "V",
}

# Optional common ambiguous codes if you want to allow them later
DNA_AMBIGUOUS: Set[str] = {"N"}  # keep minimal unless assignment says otherwise
PROTEIN_AMBIGUOUS: Set[str] = {"X"}  # unknown amino acid


@dataclass(frozen=True)
class Alphabet:
    name: str
    symbols: Set[str]

    def normalize(self, seq: str) -> str:
        return seq.strip().upper()

    def validate(self, seq: str, allow_ambiguous: bool = False) -> bool:
        s = self.normalize(seq)
        if not s:
            return False
        allowed = set(self.symbols)
        if allow_ambiguous:
            if self.name == "DNA":
                allowed |= DNA_AMBIGUOUS
            elif self.name == "PROTEIN":
                allowed |= PROTEIN_AMBIGUOUS
        return all(ch in allowed for ch in s)


DNA = Alphabet("DNA", DNA_ALPHABET)
PROTEIN = Alphabet("PROTEIN", PROTEIN_ALPHABET)


def detect_alphabet(seq: str, allow_ambiguous: bool = False) -> AlphabetType:
    """
    Guess whether a sequence is DNA or protein.

    Heuristic:
    - If sequence contains protein-only characters (not in DNA alphabet):
      * If it also has a significant proportion of A, C, G, T (DNA-like), it's mixed -> UNKNOWN
      * Otherwise, if valid protein -> PROTEIN
    - If sequence only contains characters in DNA alphabet (A, C, G, T):
      * If valid DNA -> DNA
    - Otherwise -> UNKNOWN

    Note: Since A, C, G, T are in both alphabets, sequences with only these
    characters will be classified as DNA. Sequences with protein-specific
    amino acids will be classified as PROTEIN unless they look mixed.
    """
    s = seq.strip().upper()
    if not s:
        return "UNKNOWN"

    # Characters that are in protein alphabet but not in DNA alphabet
    protein_only_chars = PROTEIN_ALPHABET - DNA_ALPHABET
    
    # Check what characters are present
    seq_chars = set(s)
    has_protein_only = bool(seq_chars & protein_only_chars)
    
    # Check validation
    is_valid_dna = DNA.validate(s, allow_ambiguous=allow_ambiguous)
    is_valid_protein = PROTEIN.validate(s, allow_ambiguous=allow_ambiguous)
    
    # If sequence contains protein-only characters
    if has_protein_only:
        # Count how many characters are DNA-only (A, C, G, T) vs protein-only
        dna_char_count = sum(1 for ch in s if ch in DNA_ALPHABET)
        protein_only_count = sum(1 for ch in s if ch in protein_only_chars)
        
        # If it has a significant mix of both types, it's ambiguous
        # (e.g., "ACGTMK" has 4 DNA chars and 2 protein-only chars)
        if dna_char_count > 0 and protein_only_count > 0:
            # Check if it looks like a short mixed sequence (ambiguous)
            # vs a longer protein sequence (where T is just threonine)
            if len(s) <= 10 and dna_char_count >= 2:
                # Short sequence with mix -> likely mixed/unknown
                return "UNKNOWN"
            # Longer sequences with T are likely protein (T = threonine)
            if is_valid_protein:
                return "PROTEIN"
            return "UNKNOWN"
        
        # Only protein-only chars (shouldn't happen in practice)
        if is_valid_protein:
            return "PROTEIN"
        return "UNKNOWN"
    
    # Sequence only contains A, C, G, T (all valid in both alphabets)
    # Default to DNA for sequences with only these characters
    if is_valid_dna:
        return "DNA"
    
    # Not valid DNA either
    return "UNKNOWN"


def validate_sequence(seq: str, alphabet: Alphabet, allow_ambiguous: bool = False) -> None:
    """Raise ValueError if sequence violates the alphabet."""
    if not alphabet.validate(seq, allow_ambiguous=allow_ambiguous):
        bad = sorted(set(seq.upper()) - alphabet.symbols)
        raise ValueError(
            f"Invalid symbols for {alphabet.name}: {bad}. "
            f"Allowed: {sorted(alphabet.symbols)}"
        )


def detect_kind(seq: str) -> str:
    """
    Detect if a sequence is DNA or protein, returning lowercase string.

    This is a convenience wrapper around detect_alphabet that returns
    lowercase "dna" or "protein" (or "unknown") for compatibility with main.py.

    Parameters
    ----------
    seq : str
        The sequence to classify.

    Returns
    -------
    str
        One of "dna", "protein", or "unknown" (all lowercase).
    """
    result = detect_alphabet(seq)
    return result.lower() if result != "UNKNOWN" else "unknown"