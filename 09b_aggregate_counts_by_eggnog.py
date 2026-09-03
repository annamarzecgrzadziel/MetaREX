#!/usr/bin/env python3

import os
from pathlib import Path
import pandas as pd
from collections import defaultdict

def ask_path(prompt, must_be_file=False, must_be_dir=False):
    path = Path(input(prompt).strip())

    if not path.exists():
        raise FileNotFoundError(f" Not found: {path}")

    if must_be_file and not path.is_file():
        raise FileNotFoundError(f" This is not a file: {path}")

    if must_be_dir and not path.is_dir():
        raise NotADirectoryError(f" This is not a directory: {path}")

    return path


def find_quant_files(base_dir: Path):
    quants = sorted(base_dir.rglob("quant.sf"))

    if not quants:
        raise RuntimeError(f" No quant.sf files found in: {base_dir}")

    return quants


def infer_sample_name(qf: Path):
    if qf.parent.name == "quant":
        return qf.parent.parent.name
    return qf.parent.name


def load_salmon_counts(quant_files):
    all_counts = {}
    sample_names = set()

    for qf in quant_files:
        sample = infer_sample_name(qf)

        if sample in sample_names:
            raise RuntimeError(f" Duplicate sample name: {sample}")

        sample_names.add(sample)

        df = pd.read_csv(qf, sep="\t", usecols=["Name", "NumReads"])
        df = df.rename(columns={"Name": "CDS", "NumReads": sample})
        all_counts[sample] = df.set_index("CDS")[sample]

    counts_df = pd.concat(all_counts.values(), axis=1)
    counts_df.columns = list(all_counts.keys())
    counts_df = counts_df.fillna(0)

    counts_df = counts_df.round().astype(int)

    return counts_df


def load_eggnog(anno_file: Path):
    header = None

    with open(anno_file) as f:
        for line in f:
            if line.startswith("#query"):
                header = line.lstrip("#").strip().split("\t")
                break

    if header is None:
        raise RuntimeError(" Header #query not found in eggNOG file")

    df = pd.read_csv(
        anno_file,
        sep="\t",
        comment="#",
        names=header,
        low_memory=False
    )

    if "query" not in df.columns:
        raise RuntimeError(" Column 'query' is still missing after parsing")

    df = df.rename(columns={"query": "CDS"})
    return df


def aggregate(counts, anno, colname):
    if colname not in anno.columns:
        raise RuntimeError(f" Column '{colname}' not found in eggNOG file")

    mapping = defaultdict(list)

    for _, row in anno.iterrows():
        cds = row["CDS"]

        if pd.isna(row[colname]):
            continue

        for term in str(row[colname]).split(","):
            term = term.strip()
            if term and term != "-":
                mapping[cds].append(term)

    agg = defaultdict(lambda: [0] * counts.shape[1])

    for cds, values in counts.iterrows():
        if cds not in mapping:
            continue

        values_list = values.tolist()

        for term in mapping[cds]:
            agg[term] = [
                agg[term][i] + int(values_list[i])
                for i in range(len(values_list))
            ]

    agg_df = pd.DataFrame.from_dict(
        agg,
        orient="index",
        columns=counts.columns
    )

    agg_df.index.name = colname
    return agg_df.sort_index()


def main():
    
    salmon_dir = ask_path(
        " Enter the directory containing Salmon results (with sample folders): ",
        must_be_dir=True
    )

    eggnog_file = ask_path(
        " Enter the eggNOG file (*.emapper.annotations): ",
        must_be_file=True
    )

    out_dir = Path(input(" Enter the output directory: ").strip())
    out_dir.mkdir(parents=True, exist_ok=True)

    quant_files = find_quant_files(salmon_dir)
    print(f" Found {len(quant_files)} quant.sf files")

    counts = load_salmon_counts(quant_files)
    print(f" {counts.shape[0]} CDS × {counts.shape[1]} samples")

    counts_out = out_dir / "09b_COUNTS_CDS_matrix.tsv"
    counts.to_csv(counts_out, sep="\t")
    print(f" Saved CDS COUNTS matrix: {counts_out}")

    eggnog = load_eggnog(eggnog_file)
    print(f" {eggnog.shape[0]} eggNOG records")

    for label, col in {
        "KO": "KEGG_ko",
        "EC": "EC",
        "PFAM": "PFAMs"
    }.items():
        print(f"\n Aggregating COUNTS → {label}")
        mat = aggregate(counts, eggnog, col)
        out_file = out_dir / f"09b_COUNTS_{label}_matrix.tsv"
        mat.to_csv(out_file, sep="\t")
        print(f" {label}: {mat.shape[0]} functions → {out_file}")

print("\n DONE")


if __name__ == "__main__":
    main()
