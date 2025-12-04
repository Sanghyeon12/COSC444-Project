# src/eval/baseline.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Any, Sequence
import csv
import random

from .randomize import randomize_sequences
from . import metrics as m


MatchList = Sequence[Mapping[str, Any]]
SequenceMap = Mapping[str, str]


@dataclass
class BaselineResult:
    """Container for baseline statistics."""
    real_counts: Dict[str, int]
    random_counts: Dict[str, List[int]]
    comparison: Dict[str, m.ComparisonMetrics]


def run_baseline(
    sequences: SequenceMap,
    scan_func: Callable[[SequenceMap], MatchList],
    trials: int = 10,
    seed: int | None = None,
) -> BaselineResult:
    """
    Run baseline evaluation.

    Args:
        sequences: Mapping from seq_id -> sequence string.
        scan_func: Function that runs the scanner and returns
                   a list of match dicts with a 'motif' key.
        trials: Number of randomization trials.
        seed: Optional RNG seed for reproducibility.

    Returns:
        BaselineResult with real counts, random counts, and comparison metrics.
    """
    rng = random.Random(seed)

    # Real data scan
    real_matches = scan_func(sequences)
    real_counts = m.count_matches_by_motif(real_matches)

    # Randomized trials
    random_counts: Dict[str, List[int]] = {motif: [] for motif in real_counts}

    for _ in range(trials):
        randomized = randomize_sequences(sequences, rng)
        trial_matches = scan_func(randomized)
        trial_counts = m.count_matches_by_motif(trial_matches)

        # make sure we track all motifs we have seen so far
        all_motifs = set(random_counts.keys()) | set(trial_counts.keys())
        for motif in all_motifs:
            if motif not in random_counts:
                random_counts[motif] = []

        # record trial counts, zero if motif not seen in this trial
        for motif in random_counts.keys():
            random_counts[motif].append(trial_counts.get(motif, 0))

    comparison = m.build_comparison(real_counts, random_counts)

    return BaselineResult(
        real_counts=dict(real_counts),
        random_counts=random_counts,
        comparison=comparison,
    )


def save_baseline_trials_csv(
    path: str | Path,
    random_counts: Mapping[str, Sequence[int]],
) -> None:
    """
    Write baseline trial counts to CSV.

    Columns: motif, trial, count
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)

    with p.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["motif", "trial", "count"])

        for motif, counts in sorted(random_counts.items()):
            for trial_idx, count in enumerate(counts, start=1):
                writer.writerow([motif, trial_idx, count])


def save_comparison_csv(
    path: str | Path,
    comparison: Mapping[str, m.ComparisonMetrics],
) -> None:
    """
    Write comparison metrics to CSV.

    Columns:
      motif, real_count, mean_random, std_random, min_random,
      max_random, lift
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)

    with p.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "motif",
                "real_count",
                "mean_random",
                "std_random",
                "min_random",
                "max_random",
                "lift",
            ]
        )

        for motif, cm in sorted(comparison.items()):
            writer.writerow(
                [
                    cm.motif,
                    cm.real_count,
                    f"{cm.mean_random:.6f}",
                    f"{cm.std_random:.6f}",
                    cm.min_random,
                    cm.max_random,
                    "inf" if cm.lift is None or cm.lift == float("inf") else f"{cm.lift:.6f}",
                ]
            )
