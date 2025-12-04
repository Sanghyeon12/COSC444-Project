"""
src/scan/report.py – Sadman (Person 3)

Legacy reporting functions (not currently used by main.py).

Note: The main.py module now contains its own reporting functions
(write_matches_txt, write_summary_txt, write_baseline_txt) that are
integrated with the current data structures. This module is kept for
reference but is not imported or used.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple, Any


def write_matches(out_path: str | Path, results: Dict[str, List[Tuple[str, int, int, str, str]]]) -> None:
    """
    Write match results to a TSV file (legacy format).

    Note: This function is not currently used. See main.py for the active reporting.
    """
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    
    with path.open("w", encoding="utf-8") as f:
        for seq_id, hits in results.items():
            for h in hits:
                line = f"{seq_id}\t{h[0]}\t{h[1]}\t{h[2]}\t{h[3]}\n"
                f.write(line)


def write_summary(out_path: str | Path, results: Dict[str, List[Any]]) -> None:
    """
    Write summary counts to a TSV file (legacy format).

    Note: This function is not currently used. See main.py for the active reporting.
    """
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    
    with path.open("w", encoding="utf-8") as f:
        for seq_id, hits in results.items():
            f.write(f"{seq_id}\t{len(hits)}\n")


def write_baseline(out_path: str | Path, results: Dict[str, List[Any]]) -> None:
    """
    Write baseline results to a TSV file (legacy format).

    Note: This function is not currently used. See main.py for the active reporting.
    """
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    
    with path.open("w", encoding="utf-8") as f:
        for seq_id, hits in results.items():
            f.write(f"{seq_id}\t{len(hits)}\n")
