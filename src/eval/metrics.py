# src/eval/metrics.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Mapping, Sequence, Any
import math
import statistics


@dataclass
class ComparisonMetrics:
    """Summary comparison for one motif."""
    motif: str
    real_count: int
    mean_random: float
    std_random: float
    min_random: int
    max_random: int
    lift: float | None


def count_matches_by_motif(
    matches: Sequence[Mapping[str, Any]],
) -> Dict[str, int]:
    """
    Count how many matches are observed for each motif.

    Each match is expected to have a 'motif' key.
    """
    counts: Dict[str, int] = {}
    for match in matches:
        motif = str(match["motif"])
        counts[motif] = counts.get(motif, 0) + 1
    return counts


def summarize_trial_counts(trial_counts: Sequence[int]) -> tuple[float, float, int, int]:
    """
    Compute mean, std, min, max for a list of counts.

    Uses population standard deviation (pstdev).
    """
    if not trial_counts:
        return 0.0, 0.0, 0, 0

    mean_val = statistics.mean(trial_counts)
    if len(trial_counts) > 1:
        std_val = statistics.pstdev(trial_counts)
    else:
        std_val = 0.0

    return mean_val, std_val, min(trial_counts), max(trial_counts)


def compute_lift(real_count: int, mean_random: float) -> float | None:
    """
    Compute enrichment (lift) = real_count / mean_random.

    Returns:
      * None if there is no sensible value (real=0 and mean_random=0)
      * +inf if mean_random is 0 but real_count > 0
    """
    if mean_random == 0.0:
        if real_count == 0:
            return None
        return math.inf

    return real_count / mean_random


def build_comparison(
    real_counts: Mapping[str, int],
    random_counts: Mapping[str, Sequence[int]],
) -> Dict[str, ComparisonMetrics]:
    """
    Build comparison metrics for each motif from real and random counts.
    """
    comparison: Dict[str, ComparisonMetrics] = {}

    motif_names = set(real_counts.keys()) | set(random_counts.keys())
    for motif in sorted(motif_names):
        real = real_counts.get(motif, 0)
        trials = list(random_counts.get(motif, []))

        mean_r, std_r, min_r, max_r = summarize_trial_counts(trials)
        lift = compute_lift(real, mean_r)

        comparison[motif] = ComparisonMetrics(
            motif=motif,
            real_count=real,
            mean_random=mean_r,
            std_random=std_r,
            min_random=min_r,
            max_random=max_r,
            lift=lift,
        )

    return comparison
