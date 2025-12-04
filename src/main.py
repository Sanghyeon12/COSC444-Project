# main.py – Sadman (Person 3)
import argparse
from io import fasta_reader
from motifs import motif_loader
from scan import scanner, report

def load_data(fasta_path, motif_path):
    sequences = fasta_reader.read_fasta(fasta_path)
    motifs = motif_loader.load_motifs(motif_path)
    return sequences, motifs

def run_scan(sequences, motifs, allow_overlap, revcomp):
    results = {}
    for seq_id, seq in sequences.items():
        hits = scanner.scan_sequence(seq_id, seq, motifs, allow_overlap, revcomp)
        results[seq_id] = hits
    return results

def main():
    p = argparse.ArgumentParser(description="Motif scanning tool")
    p.add_argument("--fasta", required=True, help="Input FASTA file")
    p.add_argument("--motifs", required=True, help="Motif JSON file")
    p.add_argument("--out", required=True, help="Output file")
    p.add_argument("--summary", action="store_true", help="Write summary instead of detailed matches")
    p.add_argument("--no-overlap", action="store_true", help="Disallow overlapping matches")
    p.add_argument("--revcomp", action="store_true", help="Scan reverse complement also")

    args = p.parse_args()

    sequences, motifs = load_data(args.fasta, args.motifs)

    results = run_scan(
        sequences,
        motifs,
        allow_overlap=not args.no_overlap,
        revcomp=args.revcomp
    )

    if args.summary:
        report.write_summary(args.out, results)
    else:
        report.write_matches(args.out, results)

    print("Done.")

if __name__ == "__main__":
    main()
"""
src/main.py

Creates 3 output files:
  - matches.txt  (full detail)
  - summary.txt  (counts/totals)
  - baseline.txt (real vs randomized control)
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Sequence, Tuple, Any
import re

from src.motifs import load_motifs
from src.io.alphabet import detect_kind
from src.io.fasta_reader import read_fasta
from src.scan.scanner import scan_sequence
from src.eval import run_baseline


@dataclass(frozen=True)
class Record:
    seq_id: str
    seq: str
    kind: str  # "dna" or "protein"

    @property
    def length(self) -> int:
        return len(self.seq)


def group_matches_by_seq(matches: Sequence[Mapping[str, Any]]) -> Dict[str, List[Mapping[str, Any]]]:
    by: Dict[str, List[Mapping[str, Any]]] = {}
    for m in matches:
        by.setdefault(str(m["seq_id"]), []).append(m)
    return by


def build_records(fasta_records: Sequence[Tuple[str, str]]) -> List[Record]:
    out: List[Record] = []
    for seq_id, seq in fasta_records:
        out.append(Record(seq_id=seq_id, seq=seq, kind=detect_kind(seq)))
    return out


def run_scan(records: Sequence[Record], dna_motifs: Mapping[str, re.Pattern], protein_motifs: Mapping[str, re.Pattern]) -> List[Dict[str, Any]]:
    matches: List[Dict[str, Any]] = []
    for rec in records:
        motif_map = dna_motifs if rec.kind == "dna" else protein_motifs
        for m in scan_sequence(rec.seq_id, rec.seq, motif_map):
            m["kind"] = rec.kind
            matches.append(m)
    return matches


def write_matches_txt(
    path: Path,
    *,
    fasta_path: Path,
    mode: str,
    dna_motifs: Mapping[str, re.Pattern],
    protein_motifs: Mapping[str, re.Pattern],
    records: Sequence[Record],
    matches: Sequence[Mapping[str, Any]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    by_seq = group_matches_by_seq(matches)

    with path.open("w", encoding="utf-8") as f:
        f.write("DNA/Protein Motif Scanner - Match Report\n")
        f.write(f"Input FASTA: {fasta_path.as_posix()}\n")
        f.write("Motif set : motifs/dna.json + motifs/protein.json\n")
        f.write(f"Mode      : {mode}\n")
        f.write("Indexing  : 1-based positions, end inclusive\n\n")

        for rec in records:
            f.write("=" * 60 + "\n")
            f.write(f"> {rec.seq_id}   kind={rec.kind.upper()}   length={rec.length}\n")
            f.write("-" * 60 + "\n")

            motif_map = dna_motifs if rec.kind == "dna" else protein_motifs
            rec_matches = by_seq.get(rec.seq_id, [])
            grouped: Dict[str, List[Mapping[str, Any]]] = {}
            for m in rec_matches:
                grouped.setdefault(str(m["motif"]), []).append(m)

            for motif_name in motif_map.keys():
                pat = motif_map[motif_name].pattern
                hits = grouped.get(motif_name, [])
                f.write(f"[{motif_name}] regex={pat}   hits={len(hits)}\n")
                if not hits:
                    f.write("  (no matches)\n\n")
                    continue
                for i, hit in enumerate(hits, 1):
                    f.write(f"  - hit#{i}  start={hit['start']}  end={hit['end']}  match={hit['match']}\n")
                f.write("\n")

        f.write("=" * 24 + " END OF MATCHES " + "=" * 24 + "\n")


def write_summary_txt(
    path: Path,
    *,
    fasta_path: Path,
    mode: str,
    dna_motifs: Mapping[str, re.Pattern],
    protein_motifs: Mapping[str, re.Pattern],
    records: Sequence[Record],
    matches: Sequence[Mapping[str, Any]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    dna_records = [r for r in records if r.kind == "dna"]
    prot_records = [r for r in records if r.kind == "protein"]

    dna_matches = [m for m in matches if m["kind"] == "dna"]
    prot_matches = [m for m in matches if m["kind"] == "protein"]

    motif_totals: Dict[str, int] = {}
    for m in matches:
        motif_totals[str(m["motif"])] = motif_totals.get(str(m["motif"]), 0) + 1

    per_seq: Dict[str, Dict[str, int]] = {r.seq_id: {} for r in records}
    for m in matches:
        sid = str(m["seq_id"])
        motif = str(m["motif"])
        per_seq[sid][motif] = per_seq[sid].get(motif, 0) + 1

    with path.open("w", encoding="utf-8") as f:
        f.write("DNA/Protein Motif Scanner - Summary\n")
        f.write(f"Input FASTA: {fasta_path.as_posix()}\n")
        f.write(f"Mode      : {mode}\n")
        f.write("Indexing  : 1-based positions, end inclusive\n\n")

        f.write("=" * 60 + "\n")
        f.write("TOTALS BY KIND\n")
        f.write("-" * 60 + "\n")
        f.write(f"DNA sequences      : {len(dna_records)}\n")
        f.write(f"Protein sequences  : {len(prot_records)}\n")
        f.write(f"Total sequences    : {len(records)}\n\n")
        f.write(f"Total matches (DNA)     : {len(dna_matches)}\n")
        f.write(f"Total matches (Protein) : {len(prot_matches)}\n")
        f.write(f"Total matches (All)     : {len(matches)}\n\n")

        f.write("=" * 60 + "\n")
        f.write("TOTALS BY MOTIF (ALL SEQUENCES)\n")
        f.write("-" * 60 + "\n")
        f.write("DNA motifs\n")
        for motif in dna_motifs.keys():
            f.write(f"  - {motif:<22} : {motif_totals.get(motif, 0)}\n")
        f.write("\nProtein motifs\n")
        for motif in protein_motifs.keys():
            f.write(f"  - {motif:<22} : {motif_totals.get(motif, 0)}\n")
        f.write("\n")

        f.write("=" * 60 + "\n")
        f.write("COUNTS PER SEQUENCE\n")
        f.write("-" * 60 + "\n")
        for rec in records:
            motif_map = dna_motifs if rec.kind == "dna" else protein_motifs
            total = 0
            parts = []
            for motif in motif_map.keys():
                c = per_seq[rec.seq_id].get(motif, 0)
                parts.append(f"{motif}={c}")
                total += c
            f.write(f"{rec.seq_id}\n  " + ", ".join(parts) + f"   | {rec.kind.upper()}_total={total}\n\n")

        f.write("=" * 24 + " END OF SUMMARY " + "=" * 24 + "\n")


def write_baseline_txt(
    path: Path,
    *,
    fasta_path: Path,
    mode: str,
    trials: int,
    seed: int,
    baseline_result,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8") as f:
        f.write("DNA/Protein Motif Scanner - Randomized Baseline Comparison\n")
        f.write(f"Input FASTA       : {fasta_path.as_posix()}\n")
        f.write("Randomization     : per-sequence shuffle (preserve composition), same lengths\n")
        f.write(f"Trials            : {trials}\n")
        f.write(f"Random seed       : {seed}\n")
        f.write(f"Mode              : {mode}\n")
        f.write("Metric            : count of motif matches\n\n")

        f.write("=" * 60 + "\n")
        f.write("MOTIF BASELINE TABLE (Real vs Randomized)\n")
        f.write("-" * 60 + "\n")
        f.write("Columns:\n")
        f.write("  real_count | random_mean | random_std | random_min | random_max | lift(real/mean)\n\n")

        for motif, cm in sorted(baseline_result.comparison.items()):
            lift_str = "inf" if cm.lift is None or cm.lift == float("inf") else f"{cm.lift:.2f}"
            f.write(
                f"{cm.motif}\n"
                f"  real_count={cm.real_count}  "
                f"random_mean={cm.mean_random:.2f}  "
                f"random_std={cm.std_random:.2f}  "
                f"min={cm.min_random}  "
                f"max={cm.max_random}  "
                f"lift={lift_str}\n"
            )

        f.write("\n" + "=" * 60 + "\n")
        f.write("INTERPRETATION (short)\n")
        f.write("-" * 60 + "\n")
        f.write(
            "If lift > 1, the motif appears more often in real sequences than in shuffled controls.\n"
            "Very large lift usually means the motif is highly structured (not just random chance).\n"
        )
        f.write("\n" + "=" * 24 + " END OF BASELINE " + "=" * 24 + "\n")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="DNA/Protein motif scanner (regex + FASTA)")
    parser.add_argument("--fasta", default="data/demo_mixed.fasta", help="Input FASTA file path")
    parser.add_argument("--outdir", default="out", help="Output directory")
    parser.add_argument("--mode", choices=["multi", "single"], default="multi", help="Scan mode")
    parser.add_argument("--motif", default=None, help="Motif name (required for single mode)")
    parser.add_argument("--trials", type=int, default=100, help="Random baseline trials")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for baseline")
    args = parser.parse_args(list(argv) if argv is not None else None)

    fasta_path = Path(args.fasta)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    dna_motifs = load_motifs("dna", compiled=True)
    protein_motifs = load_motifs("protein", compiled=True)

    if args.mode == "single":
        if not args.motif:
            raise SystemExit("--motif is required when --mode single")
        if args.motif in dna_motifs:
            dna_motifs = {args.motif: dna_motifs[args.motif]}
            protein_motifs = {}
        elif args.motif in protein_motifs:
            protein_motifs = {args.motif: protein_motifs[args.motif]}
            dna_motifs = {}
        else:
            raise SystemExit(f"Unknown motif: {args.motif!r}")

    fasta_records = read_fasta(fasta_path)
    records = build_records(fasta_records)
    matches = run_scan(records, dna_motifs, protein_motifs)

    write_matches_txt(outdir / "matches.txt", fasta_path=fasta_path, mode=args.mode.upper(),
                      dna_motifs=dna_motifs, protein_motifs=protein_motifs, records=records, matches=matches)

    write_summary_txt(outdir / "summary.txt", fasta_path=fasta_path, mode=args.mode.upper(),
                      dna_motifs=dna_motifs, protein_motifs=protein_motifs, records=records, matches=matches)

    seq_map = {r.seq_id: r.seq for r in records}

    def scan_func(seqs: Mapping[str, str]):
        temp_records = [Record(seq_id=k, seq=v, kind=detect_kind(v)) for k, v in seqs.items()]
        return run_scan(temp_records, dna_motifs, protein_motifs)

    baseline_result = run_baseline(seq_map, scan_func, trials=args.trials, seed=args.seed)
    write_baseline_txt(outdir / "baseline.txt", fasta_path=fasta_path, mode=args.mode.upper(),
                       trials=args.trials, seed=args.seed, baseline_result=baseline_result)

    print(f"Saved outputs to: {outdir.resolve()}")
    print(f"- {outdir / 'matches.txt'}")
    print(f"- {outdir / 'summary.txt'}")
    print(f"- {outdir / 'baseline.txt'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
