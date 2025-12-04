from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, List, Optional, TextIO, Union


@dataclass(frozen=True)
class FastaRecord:
    """A single FASTA record."""
    id: str
    description: str
    sequence: str

    @property
    def header(self) -> str:
        return f">{self.id} {self.description}".rstrip()


def _iter_fasta_lines(handle: TextIO) -> Iterator[str]:
    """Yield non-empty, stripped lines from a FASTA stream."""
    for raw in handle:
        line = raw.strip()
        if not line:
            continue
        yield line


def parse_fasta(handle: TextIO) -> Iterator[FastaRecord]:
    """
    Parse a FASTA stream into FastaRecord objects.

    Rules:
    - Lines starting with '>' begin a new record.
    - Sequence lines may be multiline.
    - Sequences are uppercased and whitespace removed.
    - Raises ValueError if format is invalid.
    """
    seq_id: Optional[str] = None
    desc: str = ""
    seq_chunks: List[str] = []

    for line in _iter_fasta_lines(handle):
        if line.startswith(">"):
            # Emit previous record if present
            if seq_id is not None:
                sequence = "".join(seq_chunks).upper()
                if not sequence:
                    raise ValueError(f"Empty sequence for FASTA record '{seq_id}'.")
                yield FastaRecord(seq_id, desc, sequence)

            # Start new record
            header = line[1:].strip()
            if not header:
                raise ValueError("FASTA header line '>' must be followed by an identifier.")
            parts = header.split(maxsplit=1)
            seq_id = parts[0]
            desc = parts[1] if len(parts) > 1 else ""
            seq_chunks = []
        else:
            if seq_id is None:
                raise ValueError("Found sequence data before first FASTA header ('>').")
            seq_chunks.append(line.replace(" ", "").replace("\t", ""))

    # Emit last record
    if seq_id is not None:
        sequence = "".join(seq_chunks).upper()
        if not sequence:
            raise ValueError(f"Empty sequence for FASTA record '{seq_id}'.")
        yield FastaRecord(seq_id, desc, sequence)


def read_fasta(
    source: Union[str, Path, TextIO]
) -> List[FastaRecord]:
    """
    Read FASTA records from a file path or an already-open handle.

    Returns a list of FastaRecord objects.
    For tuple format (seq_id, sequence), use: [(r.id, r.sequence) for r in read_fasta(...)]
    """
    if hasattr(source, "read"):
        return list(parse_fasta(source))  # type: ignore[arg-type]

    path = Path(source)
    with path.open("r", encoding="utf-8") as f:
        return list(parse_fasta(f))
