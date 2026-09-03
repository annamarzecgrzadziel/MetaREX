#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import re
import pandas as pd

def ask_file_or_dir(prompt, pattern):
    while True:
        path = input(prompt).strip().strip('"').strip("'")

        if not path:
            print(" Empty path – try again")
            continue

        if not os.path.exists(path):
            print(f" Path does not exist: {path}")
            continue

        if os.path.isfile(path):
            return path

        if os.path.isdir(path):
            matches = sorted([
                os.path.join(path, f)
                for f in os.listdir(path)
                if re.search(pattern, f)
            ])

            if not matches:
                print(f" No files matching pattern: {pattern}")
                continue

            if len(matches) == 1:
                print(f" File found: {matches[0]}")
                return matches[0]

            print(" Multiple files found:")
            for i, f in enumerate(matches, 1):
                print(f"  [{i}] {os.path.basename(f)}")

            try:
                idx = int(input(" Select file number: "))
                return matches[idx - 1]
            except Exception:
                print(" Invalid selection")
                continue


def split_features(val):
    if pd.isna(val):
        return []
    val = str(val).strip()
    if val in ("", "-", "NA", "None"):
        return []
    parts = re.split(r"[,\s;|]+", val)
    return [p for p in parts if p and p != "-"]


def load_tpm(tpm_file):
    df = pd.read_csv(tpm_file, sep="\t")
    df.rename(columns={df.columns[0]: "GeneID"}, inplace=True)
    sample_cols = list(df.columns[1:])
    df[sample_cols] = df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    return df, sample_cols


def load_eggnog(anno_file):
    
    import tempfile
    import sys

    cleaned_lines = []
    header_fixed = False

    with open(anno_file) as f:
        for line in f:
            if line.startswith("#query"):
                cleaned_lines.append(line.lstrip("#"))
                header_fixed = True
            elif line.startswith("##"):
                continue
            else:
                cleaned_lines.append(line)

    if not header_fixed:
        sys.exit(" Header '#query' not found in eggNOG file")

    with tempfile.NamedTemporaryFile("w", delete=False) as tmp:
        tmp.writelines(cleaned_lines)
        tmp_path = tmp.name

    df = pd.read_csv(tmp_path, sep="\t", low_memory=False)

    if "query" not in df.columns:
        sys.exit(
            " Column 'query' is still missing.\n"
            f"Available columns:\n{df.columns.tolist()}"
        )

    print(" Detected CDS ID column: query")

    out = pd.DataFrame()
    out["GeneID"] = df["query"].astype(str)

    out["KO"] = df["KEGG_ko"] if "KEGG_ko" in df.columns else ""
    out["EC"] = df["EC"] if "EC" in df.columns else ""
    out["PFAM"] = df["PFAMs"] if "PFAMs" in df.columns else ""

    return out




def explode_mapping(map_df, colname):
    rows = []
    for gid, feats in zip(map_df["GeneID"], map_df[colname]):
        for f in split_features(feats):
            rows.append((gid, f))
    return pd.DataFrame(rows, columns=["GeneID", "Feature"]).drop_duplicates()


def aggregate_tpm(tpm_df, sample_cols, mapping_df):
    merged = mapping_df.merge(tpm_df, on="GeneID", how="inner")
    if merged.empty:
        return pd.DataFrame(columns=["Feature"] + sample_cols)
    agg = merged.groupby("Feature", as_index=False)[sample_cols].sum()
    agg["__sum"] = agg[sample_cols].sum(axis=1)
    agg = agg.sort_values("__sum", ascending=False).drop(columns="__sum")
    return agg


def main():
    print("\n metaTP – Step 09: TPM aggregation (CDS → KO / EC / Pfam)")
    print("========================================================\n")

    tpm_file = ask_file_or_dir(
        " Enter TPM file OR directory (TPM1.0_samples1.tsv): ",
        pattern=r"TPM.*\.tsv$"
    )

    eggnog_file = ask_file_or_dir(
        " Enter eggNOG annotations file OR directory (*.emapper.annotations): ",
        pattern=r"\.emapper\.annotations$"
    )

    out_dir = input(" Enter output directory: ").strip()
    os.makedirs(out_dir, exist_ok=True)

    tpm_df, sample_cols = load_tpm(tpm_file)
    print(f" {tpm_df.shape[0]} CDS, {len(sample_cols)} samples")

    eggnog_df = load_eggnog(eggnog_file)
    print(f" {eggnog_df.shape[0]} eggNOG records")

    overlap = len(set(tpm_df["GeneID"]) & set(eggnog_df["GeneID"]))
    print(f" Overlap CDS (TPM ∩ eggNOG): {overlap}")

    ko_map = explode_mapping(eggnog_df, "KO")
    ec_map = explode_mapping(eggnog_df, "EC")
    pf_map = explode_mapping(eggnog_df, "PFAM")

    print(f"  KO:   {ko_map.shape[0]}")
    print(f"  EC:   {ec_map.shape[0]}")
    print(f"  PFAM: {pf_map.shape[0]}")

    ko_mat = aggregate_tpm(tpm_df, sample_cols, ko_map)
    ec_mat = aggregate_tpm(tpm_df, sample_cols, ec_map)
    pf_mat = aggregate_tpm(tpm_df, sample_cols, pf_map)

    ko_out = os.path.join(out_dir, "09_TPM_KO_matrix.tsv")
    ec_out = os.path.join(out_dir, "09_TPM_EC_matrix.tsv")
    pf_out = os.path.join(out_dir, "09_TPM_PFAM_matrix.tsv")

    ko_mat.to_csv(ko_out, sep="\t", index=False)
    ec_mat.to_csv(ec_out, sep="\t", index=False)
    pf_mat.to_csv(pf_out, sep="\t", index=False)

    print("\n Files saved:")
    print(f"   KO:   {ko_out}   ({ko_mat.shape[0]} functions)")
    print(f"    EC:   {ec_out}   ({ec_mat.shape[0]} functions)")
    print(f"   PFAM: {pf_out}   ({pf_mat.shape[0]} functions)")

    print("\n Step 09 completed")
    print(" Next step: 10_ALDEx2_metaTP_full.R\n")


if __name__ == "__main__":
    main()

