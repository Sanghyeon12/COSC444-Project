from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Dict, Mapping

from .iupac import iupac_to_regex


def _motifs_dir() -> Path:
    """Return the path to the top-level ``motifs/`` directory."""
    here = Path(__file__).resolve()
    for parent in here.parents:
        candidate = parent / "motifs"
        if candidate.is_dir() and (candidate / "dna.json").is_file() and (candidate / "protein.json").is_file():
            return candidate
    raise RuntimeError("Could not locate top-level 'motifs' directory containing dna.json / protein.json")


def _load_raw_json(kind: str) -> Dict[str, str]:
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
