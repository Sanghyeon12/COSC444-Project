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
    - If all chars are valid DNA -> DNA
    - Else if all chars are valid protein -> PROTEIN
    - Else UNKNOWN

    This works well for mixed FASTA datasets.
    """
    s = seq.strip().upper()
    if not s:
        return "UNKNOWN"

    if DNA.validate(s, allow_ambiguous=allow_ambiguous):
        return "DNA"
    if PROTEIN.validate(s, allow_ambiguous=allow_ambiguous):
        return "PROTEIN"
    return "UNKNOWN"


def validate_sequence(seq: str, alphabet: Alphabet, allow_ambiguous: bool = False) -> None:
    """Raise ValueError if sequence violates the alphabet."""
    if not alphabet.validate(seq, allow_ambiguous=allow_ambiguous):
        bad = sorted(set(seq.upper()) - alphabet.symbols)
        raise ValueError(
            f"Invalid symbols for {alphabet.name}: {bad}. "
            f"Allowed: {sorted(alphabet.symbols)}"
        )
