"""
src/scan/scanner.py

Scanner module for finding motif matches in DNA and protein sequences.
Supports both single-motif and multi-motif scanning modes.
"""

from __future__ import annotations

import re
from typing import Dict, List, Mapping, Any


def scan_sequence(
    seq_id: str,
    sequence: str,
    motif_map: Mapping[str, re.Pattern],
) -> List[Dict[str, Any]]:
    """
    Scan a sequence for motif matches using compiled regex patterns.

    Parameters
    ----------
    seq_id : str
        Identifier for the sequence.
    sequence : str
        The DNA or protein sequence to scan.
    motif_map : Mapping[str, re.Pattern]
        Dictionary mapping motif names to compiled regex patterns.

    Returns
    -------
    List[Dict[str, Any]]
        List of match dictionaries, each containing:
        - seq_id: sequence identifier
        - motif: motif name
        - start: 0-based start position
        - end: 0-based end position (exclusive)
        - match: matched substring

    Notes
    -----
    - Uses overlapping matches (all matches are reported)
    - Positions are 0-based (will be converted to 1-based in reporting)
    - Supports multi-motif scanning (all motifs in motif_map are searched)
    """
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
