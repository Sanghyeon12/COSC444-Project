# src/eval/randomize.py
from __future__ import annotations

from typing import Dict, Mapping
import random


SequenceMap = Dict[str, str]


def randomize_sequence(seq: str, rng: random.Random | None = None) -> str:
    """
    Return a randomized version of a sequence.

    The function preserves the multiset of characters
    (same letters, different order).
    """
    if rng is None:
        rng = random.Random()

    chars = list(seq)
    rng.shuffle(chars)
    return "".join(chars)


def randomize_sequences(
    sequences: Mapping[str, str],
    rng: random.Random | None = None,
) -> SequenceMap:
    """
    Randomize all sequences in a mapping {seq_id: sequence}.

    Character composition of each sequence is preserved.
    """
    if rng is None:
        rng = random.Random()

    randomized: SequenceMap = {}
    for seq_id, seq in sequences.items():
        randomized[seq_id] = randomize_sequence(seq, rng)
    return randomized
