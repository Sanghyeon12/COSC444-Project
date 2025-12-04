"""
src.io package

FASTA file reading and alphabet detection utilities.
"""

from __future__ import annotations

from .fasta_reader import read_fasta, parse_fasta, FastaRecord
from .alphabet import (
    detect_alphabet,
    detect_kind,
    validate_sequence,
    DNA,
    PROTEIN,
    Alphabet,
    DNA_ALPHABET,
    PROTEIN_ALPHABET,
)

__all__ = [
    "read_fasta",
    "parse_fasta",
    "FastaRecord",
    "detect_alphabet",
    "detect_kind",
    "validate_sequence",
    "DNA",
    "PROTEIN",
    "Alphabet",
    "DNA_ALPHABET",
    "PROTEIN_ALPHABET",
]

