from __future__ import annotations

import re
from typing import Dict, List, Mapping, Any


def scan_sequence(
    seq_id: str,
    sequence: str,
    motif_map: Mapping[str, re.Pattern],
) -> List[Dict[str, Any]]:
    matches: List[Dict[str, Any]] = []

    for motif_name, pattern in motif_map.items():
        # Find all overlapping matches for this motif
        for match_obj in pattern.finditer(sequence):
            matches.append({
                "seq_id": seq_id,
                "motif": motif_name,
                "start": match_obj.start() + 1,  # Convert to 1-based (inclusive)
                "end": match_obj.end(),  # 1-based end position (inclusive)
                "match": match_obj.group(0),
            })

    return matches
