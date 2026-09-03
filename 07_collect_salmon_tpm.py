#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
import sys


def ask_path(prompt, must_exist=True, is_dir=True):
    path = Path(input(prompt).strip())
    if must_exist and not path.exists():
        print(f"Path does not exist: {path}")
        sys.exit(1)
    if is_dir and must_exist and not path.is_dir():
        print(f"This is not a directory: {path}")
        sys.exit(1)
    return path


def ask_float(prompt, default):
    val = input(f"{prompt} [default: {default}]: ").strip()
    return float(val) if val else default


def ask_int(prompt, default):
    val = input(f"{prompt} [default: {default}]: ").strip()
    return int(val) if val else default


def main():
    print("\nCDS FILTERING BY EXPRESSION (Salmon TPM)\n")

    input_dir = ask_path(
        "Enter the directory containing Salmon results (e.g. 06_quant_salmon):\n> "
    )

    quant_files = sorted(input_dir.glob("*/quant.sf"))
    if not quant_files:
        print("No quant.sf files found in the specified directory")
        sys.exit(1)

    print(f"Found {len(quant_files)} samples\n")

    output_dir = ask_path(
        "Enter the directory for saving results:\n> ",
        must_exist=False
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    tpm_threshold = ask_float("Minimum TPM threshold", 1.0)
    min_samples = ask_int("Minimum number of samples meeting the threshold", 1)

    tpm_all_path = output_dir / "salmon_TPM_all.tsv"
    tpm_filt_path = output_dir / f"salmon_TPM_TPM{tpm_threshold}_samples{min_samples}.tsv"

    dfs = []
    sample_names = []

    for qf in quant_files:
        sample = qf.parent.name

        if sample in sample_names:
            print(f"Duplicate sample name: {sample}")
            print(f"   Check the directory structure in: {input_dir}")
            sys.exit(1)

        sample_names.append(sample)

        try:
            df = pd.read_csv(qf, sep="\t", usecols=["Name", "TPM"])
        except Exception as e:
            print(f" Failed to read file: {qf}")
            print(f"   Details: {e}")
            sys.exit(1)

        df = df.rename(columns={"TPM": sample})
        dfs.append(df)

    merged = dfs[0]
    for df in dfs[1:]:
        merged = merged.merge(df, on="Name", how="outer")

    merged = merged.fillna(0)
    merged.to_csv(tpm_all_path, sep="\t", index=False)

    expr = merged.iloc[:, 1:]
    mask = (expr >= tpm_threshold).sum(axis=1) >= min_samples
    filtered = merged.loc[mask].copy()

    filtered.to_csv(tpm_filt_path, sep="\t", index=False)

     print("\n FILTERING COMPLETED\n")
    print(f" TPM matrix (complete): {tpm_all_path}")
    print(f" TPM matrix (filtered): {tpm_filt_path}\n")
    print(f" Samples:       {len(sample_names)}")
    print(f" CDS before:    {merged.shape[0]}")
    print(f" CDS after:     {filtered.shape[0]}")
    print(f" Removed:       {merged.shape[0] - filtered.shape[0]}")
    print("\n Ready for downstream analysis (functions / DE / visualization)\n")


if __name__ == "__main__":
    main()
