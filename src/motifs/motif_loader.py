"""
motif_loader.py - Sanghyeon (Person 2)

Utilities for loading built-in DNA / protein motifs from the JSON files in
the top-level ``motifs/`` directory.

Each JSON file should map motif names to *IUPAC / regex* pattern strings.
We translate any IUPAC codes to plain regular expressions and optionally
compile them with :mod:`re`.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Dict, Mapping

from .iupac import iupac_to_regex


def _motifs_dir() -> Path:
    """Return the path to the top-level ``motifs/`` directory.

    This walks up from the current file until it finds a ``motifs`` folder.
    It keeps the project layout flexible and works when running tests from
    the repository root.
    """
    here = Path(__file__).resolve()
    for parent in here.parents:
        candidate = parent / "motifs"
        if candidate.is_dir():
            return candidate
    # If we get here something is wrong with the checkout / install
    raise RuntimeError("Could not locate 'motifs' directory relative to motifs package")


def _load_raw_json(kind: str) -> Dict[str, str]:
    """Load the raw motif definitions from ``motifs/<kind>.json``.

    Parameters
    ----------
    kind:
        Either ``"dna"`` or ``"protein"`` (case-insensitive).

    Returns
    -------
    dict
        Mapping of motif name → pattern string (still IUPAC / regex).

    Raises
    ------
    FileNotFoundError
        If the JSON file cannot be found.
    ValueError
        If ``kind`` is not recognised.
    """
    kind_lower = kind.lower()
    if kind_lower not in {"dna", "protein"}:
        raise ValueError(f"Unknown motif kind: {kind!r} (expected 'dna' or 'protein')")

    path = _motifs_dir() / f"{kind_lower}.json"
    with path.open("r", encoding="utf-8") as fh:
        data = json.load(fh)

    # JSON is expected to be a simple name → pattern mapping
    if not isinstance(data, dict):
        raise ValueError(f"Motif JSON {path} must contain an object mapping names to patterns")
    return {str(name): str(pattern) for name, pattern in data.items()}


def load_motifs(kind: str, *, compiled: bool = True) -> Mapping[str, re.Pattern | str]:
    """Load built-in motifs for the given alphabet.

    Parameters
    ----------
    kind:
        Either ``"dna"`` or ``"protein"``.
    compiled:
        If ``True`` (default), return a mapping of motif name → compiled
        :class:`re.Pattern`. If ``False``, return motif name → regex *string*.

    Returns
    -------
    Mapping[str, Pattern | str]
        The loaded motif definitions.

    Notes
    -----
    - The JSON values are treated as IUPAC / regex strings. Any recognised
      IUPAC codes are expanded using :func:`iupac_to_regex`, but existing
      regex constructs like ``[AT]`` or ``[^P]`` are preserved.
    """
    raw = _load_raw_json(kind)
    kind_lower = kind.lower()

    # First expand any IUPAC codes
    regex_strings: Dict[str, str] = {
        name: iupac_to_regex(pattern, kind_lower)
        for name, pattern in raw.items()
    }

    if not compiled:
        return regex_strings

    compiled_dict: Dict[str, re.Pattern] = {
        name: re.compile(regex)
        for name, regex in regex_strings.items()
    }
    return compiled_dict
