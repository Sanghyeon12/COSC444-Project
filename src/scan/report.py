from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple, Any


def write_matches(out_path: str | Path, results: Dict[str, List[Tuple[str, int, int, str, str]]]) -> None:
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    
    with path.open("w", encoding="utf-8") as f:
        for seq_id, hits in results.items():
            for h in hits:
                line = f"{seq_id}\t{h[0]}\t{h[1]}\t{h[2]}\t{h[3]}\n"
                f.write(line)


def write_summary(out_path: str | Path, results: Dict[str, List[Any]]) -> None:
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    
    with path.open("w", encoding="utf-8") as f:
        for seq_id, hits in results.items():
            f.write(f"{seq_id}\t{len(hits)}\n")


def write_baseline(out_path: str | Path, results: Dict[str, List[Any]]) -> None:
    path = Path(out_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    
    with path.open("w", encoding="utf-8") as f:
        for seq_id, hits in results.items():
            f.write(f"{seq_id}\t{len(hits)}\n")
