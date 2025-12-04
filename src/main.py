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
