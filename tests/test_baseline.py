# tests/test_baseline.py
from __future__ import annotations

import math
import random

from src.eval import randomize, metrics, baseline


def test_randomize_sequence_preserves_characters():
    seq = "AABBC"
    rng = random.Random(123)

    shuffled = randomize.randomize_sequence(seq, rng)

    # Same multiset of characters
    assert sorted(shuffled) == sorted(seq)
    # With this seed, it should usually change order
    assert shuffled != seq


def test_randomize_sequences_deterministic_with_seed():
    seqs = {"s1": "AAAA", "s2": "CCGG"}

    rng1 = random.Random(42)
    rng2 = random.Random(42)

    r1 = randomize.randomize_sequences(seqs, rng1)
    r2 = randomize.randomize_sequences(seqs, rng2)

    assert r1 == r2
    assert set(r1.keys()) == set(seqs.keys())


def test_count_matches_by_motif():
    matches = [
        {"motif": "motif1"},
        {"motif": "motif1"},
        {"motif": "motif2"},
    ]
    counts = metrics.count_matches_by_motif(matches)

    assert counts["motif1"] == 2
    assert counts["motif2"] == 1


def test_build_comparison_and_lift():
    real_counts = {"M1": 10, "M2": 0}
    random_counts = {"M1": [5, 5, 5], "M2": [0, 0, 0]}

    comp = metrics.build_comparison(real_counts, random_counts)

    m1 = comp["M1"]
    assert m1.mean_random == 5
    assert math.isclose(m1.lift, 2.0)

    m2 = comp["M2"]
    assert m2.mean_random == 0
    # both zero -> lift None
    assert m2.lift is None


def _fake_scan_count_A(sequences):
    """
    Fake scanner: every 'A' is a match of motif 'A'.
    Returns list of dicts with 'motif' key like a real scanner might.
    """
    matches = []
    for seq_id, seq in sequences.items():
        for i, ch in enumerate(seq):
            if ch == "A":
                matches.append(
                    {
                        "seq_id": seq_id,
                        "motif": "A",
                        "start": i,
                        "end": i + 1,
                        "match_text": "A",
                    }
                )
    return matches


def test_run_baseline_with_fake_scanner():
    # Randomization preserves characters, so the number of 'A's
    # should be constant across trials.
    sequences = {
        "s1": "AAAA",
        "s2": "CCCC",
    }

    result = baseline.run_baseline(
        sequences,
        scan_func=_fake_scan_count_A,
        trials=5,
        seed=1,
    )

    # real: 4 A's total
    assert result.real_counts["A"] == 4

    # Each trial should also report 4 A's
    assert result.random_counts["A"] == [4, 4, 4, 4, 4]

    comp = result.comparison["A"]
    assert comp.mean_random == 4
    assert math.isclose(comp.lift, 1.0)
