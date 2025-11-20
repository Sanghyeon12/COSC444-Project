# src/eval/__init__.py
from .randomize import randomize_sequence, randomize_sequences
from .metrics import (
    ComparisonMetrics,
    count_matches_by_motif,
    summarize_trial_counts,
    compute_lift,
    build_comparison,
)
from .baseline import (
    BaselineResult,
    run_baseline,
    save_baseline_trials_csv,
    save_comparison_csv,
)

__all__ = [
    "randomize_sequence",
    "randomize_sequences",
    "ComparisonMetrics",
    "count_matches_by_motif",
    "summarize_trial_counts",
    "compute_lift",
    "build_comparison",
    "BaselineResult",
    "run_baseline",
    "save_baseline_trials_csv",
    "save_comparison_csv",
]
