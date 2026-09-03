#!/usr/bin/env python3

from Bio import SeqIO
import pandas as pd
import argparse

def get_stats(cds_file):
    lengths = [len(rec.seq) for rec in SeqIO.parse(cds_file, "fasta")]
    sample = cds_file.split("/")[-1].replace(".transcripts.fa.transdecoder.cds", "")
    return {
        "Sample": sample,
        "CDS_count": len(lengths),
        "Min_length": min(lengths),
        "Max_length": max(lengths),
        "Mean_length": round(sum(lengths) / len(lengths), 2),
        "Median_length": round(pd.Series(lengths).median(), 2),
        "Total_length": sum(lengths)
    }

def main():
    parser = argparse.ArgumentParser(description="Comparison of CDS statistics from TransDecoder")
    parser.add_argument("cds_files", nargs="+", help=".cds files to compare")
    parser.add_argument("-o", "--output", default="cds_comparison.tsv", help="Output TSV file")

    args = parser.parse_args()
    stats = [get_stats(f) for f in args.cds_files]
    df = pd.DataFrame(stats)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"Comparison saved to: {args.output}")

if __name__ == "__main__":
    main()
