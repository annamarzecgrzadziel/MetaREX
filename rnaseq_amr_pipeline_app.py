#!/usr/bin/env python3
"""
Local interactive RNA-seq AMR metatranscriptomics pipeline app.

Modes:
  python3 rnaseq_amr_pipeline_app.py serve --host 127.0.0.1 --port 8791
  python3 rnaseq_amr_pipeline_app.py run --config /path/to/pipeline_config.json

The app uses only Python standard library modules. Heavy work is delegated to
local conda environments, databases and tools already present under /data.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import json
import os
import random
import re
import shlex
import shutil
import subprocess
import sys
import threading
import time
import traceback
from collections import defaultdict
from dataclasses import dataclass
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlparse


APP_DIR = Path(__file__).resolve().parent

DEFAULT_ENVS = {
    "qc": os.environ.get("METAREX_ENV_QC", ""),
    "meta": os.environ.get("METAREX_ENV_META", ""),
    "quast": os.environ.get("METAREX_ENV_QUAST", ""),
    "eggnog": os.environ.get("METAREX_ENV_EGGNOG", ""),
    "card": os.environ.get("METAREX_ENV_CARD", ""),
    "r": os.environ.get("METAREX_ENV_R", ""),
}

DEFAULT_DBS = {
    "rrna_bowtie2": os.environ.get("METAREX_RRNA_INDEX", ""),
    "eggnog_data": os.environ.get("METAREX_EGGNOG_DATA", ""),
}

_SCRIPT_DIR = APP_DIR / "scripts"
SOURCE_SCRIPTS = {
    name: str(_SCRIPT_DIR / filename)
    for name, filename in {
        "estimate_rrna": "00_estimate_rRNA_content.sh",
        "fastp": "01_fastp_trim_QC.sh",
        "remove_rrna": "02_remove_rRNA_bowtie2.sh",
        "rrna_stats": "02a_remove_rRNA_stats.sh",
        "megahit": "03a_megahit.sh",
        "collect_contigs": "03b_contigs_all.sh",
        "quast": "03c_quast_metaquast.sh",
        "transdecoder": "04_transdecoder.sh",
        "collect_cds": "04a_transdecoder_all_cds.sh",
        "cds_stats": "05_compare_cds_stats.py",
        "salmon": "06_quant_salmon.sh",
        "salmon_tpm": "07_collect_salmon_tpm.py",
        "eggnog": "08_eggnog_mapper2.sh",
        "aggregate_tpm": "09a_aggregate_tpm_by_eggnog.py",
        "aggregate_counts": "09b_aggregate_counts_by_eggnog.py",
        "compare_ko_aldex2": "10a_DEG_KO_ALDEx2.R",
        "compare_ko_deseq2": "10b_DEG_KO_DeSeq2.R",
        "compare_pathway_aldex2": "10a_DEG_KEGG_Pathway_ALDEx2.R",
        "compare_pathway_deseq2": "10b_DEG_KEGG_Pathway_DeSeq2.R",
        "kegg_pathway": "11_KO_to_KEGG_Pathway.R",
        "rnaseq_overview": "13_RNAseq_one_group_analysis.R",
        "amr_ko": "14_AMR_KO.R",
        "card_rgi": "15_CARD_RGI_protein.sh",
        "card_analysis": "15a_CARD_analysis.R",
        "card_ko": "16_integracja_CARD_KO.R",
    }.items()
}

FASTQ_RE = re.compile(r"\.(fastq|fq)(\.gz)?$", re.IGNORECASE)
FASTA_RE = re.compile(r"\.(fa|fasta|fna)(\.gz)?$", re.IGNORECASE)
DERIVED_DIRS = {
    "fastq_clean",
    "rRNA_removed",
    "00_test_run_downsampled",
    "01_fastp_trim_QC",
    "02_rRNA_removed",
    "03_assembly_megahit",
    "06_quant_salmon",
}
DERIVED_FASTQ_RE = re.compile(r"(_clean|\.rRNAfree_|\.test)(?:[._-]|$)", re.IGNORECASE)

PIPELINE_STEPS = [
    ("estimate_rrna", "Estimate rRNA content"),
    ("trim_qc", "fastp trim + MultiQC"),
    ("remove_rrna", "Remove rRNA with Bowtie2"),
    ("rrna_stats", "rRNA removal stats"),
    ("megahit", "MEGAHIT transcript assembly"),
    ("collect_contigs", "Collect contigs"),
    ("quast", "QUAST / MetaQUAST de novo report"),
    ("transdecoder", "TransDecoder CDS/PEP"),
    ("collect_cds", "Collect CDS and proteins"),
    ("cds_stats", "CDS statistics"),
    ("salmon", "Salmon quantification"),
    ("salmon_matrix", "Salmon TPM/count matrices"),
    ("eggnog", "eggNOG-mapper"),
    ("aggregate", "Aggregate KO/EC/PFAM"),
    ("kegg_pathways", "KO to KEGG pathways"),
    ("comparative", "Comparative KO/Pathway analysis"),
    ("rnaseq_overview", "RNA-seq overview plots"),
    ("amr_ko", "AMR from KO matrix"),
    ("card_rgi", "CARD/RGI on expressed contigs"),
    ("card_summary", "CARD summary tables"),
    ("card_ko_integration", "CARD x KO integration"),
    ("multiqc", "MultiQC stage reports"),
    ("summary", "Final HTML index"),
]


@dataclass
class Sample:
    sample: str
    r1: Path | None = None
    r2: Path | None = None


def q(value: Any) -> str:
    return shlex.quote(str(value))


def mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def strip_seq_ext(name: str) -> str:
    return re.sub(r"\.(fastq|fq|fasta|fa|fna)(\.gz)?$", "", name, flags=re.I)


def open_text_maybe_gzip(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return path.open(errors="replace")


def is_derived_fastq(path: Path, input_dir: Path) -> bool:
    try:
        parts = [p.lower() for p in path.relative_to(input_dir).parts[:-1]]
    except ValueError:
        parts = [p.lower() for p in path.parts[:-1]]
    if any(part in DERIVED_DIRS for part in parts):
        return True
    return bool(DERIVED_FASTQ_RE.search(path.name))


def sample_from_r1(path: Path) -> str:
    stem = strip_seq_ext(path.name)
    patterns = [
        r"(_L\d{3})?_R1(_\d{3})?$",
        r"(_L\d{3})?_1(_\d{3})?$",
        r"[._-]R1$",
        r"[._-]1$",
    ]
    for pattern in patterns:
        cleaned = re.sub(pattern, "", stem, flags=re.I)
        if cleaned != stem:
            return cleaned
    return stem


def r2_candidates(r1: Path) -> list[Path]:
    replacements = [
        ("_R1_", "_R2_"),
        ("_R1.", "_R2."),
        ("_R1", "_R2"),
        ("_1.", "_2."),
        (".R1.", ".R2."),
        ("-R1-", "-R2-"),
    ]
    return [r1.with_name(r1.name.replace(old, new, 1)) for old, new in replacements if old in r1.name]


def discover_samples(input_dir: Path, recursive: bool = True) -> list[Sample]:
    if not input_dir.exists() or not input_dir.is_dir():
        raise RuntimeError(f"Input directory does not exist: {input_dir}")
    globber = input_dir.rglob if recursive else input_dir.glob
    raw_fastqs = sorted(p for p in globber("*") if p.is_file() and FASTQ_RE.search(p.name))
    fastqs = [p for p in raw_fastqs if not is_derived_fastq(p, input_dir)]
    if not fastqs:
        fastqs = raw_fastqs

    samples: dict[str, Sample] = {}
    used: set[Path] = set()
    for fq in fastqs:
        if fq in used:
            continue
        lower = fq.name.lower()
        if not re.search(r"([._-]r?1)([._-]|\.)", lower) and "_r1" not in lower:
            continue
        r2 = next((c for c in r2_candidates(fq) if c.exists()), None)
        if not r2:
            continue
        sample = sample_from_r1(fq)
        samples[sample] = Sample(sample=sample, r1=fq, r2=r2)
        used.add(fq)
        used.add(r2)
    return sorted(samples.values(), key=lambda item: item.sample)


def sample_table(samples: list[Sample], outdir: Path) -> Path:
    path = outdir / "00_samples.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "r1", "r2"])
        for sample in samples:
            writer.writerow([sample.sample, sample.r1 or "", sample.r2 or ""])
    return path


def load_sample_table(path: Path) -> list[Sample]:
    if not path.exists():
        return []
    samples: list[Sample] = []
    with path.open(newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            name = (row.get("sample") or "").strip()
            if name:
                samples.append(
                    Sample(
                        sample=name,
                        r1=Path(row["r1"]) if row.get("r1") else None,
                        r2=Path(row["r2"]) if row.get("r2") else None,
                    )
                )
    return samples


def fastq_record(handle: Any) -> str | None:
    lines = [handle.readline() for _ in range(4)]
    if not lines[0]:
        return None
    if any(line == "" for line in lines):
        raise RuntimeError("Truncated FASTQ record encountered during downsampling")
    return "".join(lines)


def write_random_paired_fastq_subset(r1: Path, r2: Path, out_r1: Path, out_r2: Path, target_reads: int, seed: int) -> int:
    rng = random.Random(seed)
    selected: list[tuple[str, str]] = []
    seen = 0
    with open_text_maybe_gzip(r1) as h1, open_text_maybe_gzip(r2) as h2:
        while True:
            rec1 = fastq_record(h1)
            rec2 = fastq_record(h2)
            if rec1 is None and rec2 is None:
                break
            if rec1 is None or rec2 is None:
                raise RuntimeError(f"Paired FASTQ files have different read counts: {r1} {r2}")
            seen += 1
            if len(selected) < target_reads:
                selected.append((rec1, rec2))
            else:
                idx = rng.randint(0, seen - 1)
                if idx < target_reads:
                    selected[idx] = (rec1, rec2)
    mkdir(out_r1.parent)
    with gzip.open(out_r1, "wt") as h1, gzip.open(out_r2, "wt") as h2:
        for rec1, rec2 in selected:
            h1.write(rec1)
            h2.write(rec2)
    return len(selected)


def count_fastq_reads(path: Path) -> int:
    lines = 0
    with open_text_maybe_gzip(path) as handle:
        for _line in handle:
            lines += 1
    return lines // 4


def fasta_lengths(path: Path) -> list[int]:
    lengths: list[int] = []
    current = 0
    with open_text_maybe_gzip(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
            else:
                current += len(line.strip())
    if current:
        lengths.append(current)
    return lengths


def copy_fasta_with_sample_prefix(src: Path, dst: Path, sample: str) -> None:
    prefix = f"{sample}__"
    with open_text_maybe_gzip(src) as inp, dst.open("w") as out:
        for line in inp:
            if line.startswith(">"):
                out.write(">" + prefix + line[1:])
            else:
                out.write(line)


def median(values: list[int]) -> float:
    if not values:
        return 0.0
    vals = sorted(values)
    mid = len(vals) // 2
    if len(vals) % 2:
        return float(vals[mid])
    return (vals[mid - 1] + vals[mid]) / 2.0


def bowtie2_index_exists(prefix: str) -> bool:
    base = Path(prefix)
    return any((base.parent / f"{base.name}.{suffix}").exists() for suffix in ["1.bt2", "1.bt2l"])


def parse_bowtie2_log(path: Path) -> dict[str, int]:
    single_total = 0
    single_aligned_one = 0
    single_aligned_many = 0
    single_unaligned = 0
    paired_total = 0
    paired_aligned_one = 0
    paired_aligned_many = 0
    paired_unaligned = 0
    if not path.exists():
        return {"total": 0, "aligned": 0, "unaligned": 0}
    with path.open(errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            m = re.match(r"(\d+)\s+\([^)]+\)\s+were paired; of these:", stripped)
            if m:
                paired_total = int(m.group(1))
                continue
            m = re.match(r"(\d+)\s+\([^)]+\)\s+aligned concordantly 0 times", stripped)
            if m:
                paired_unaligned = int(m.group(1))
                continue
            m = re.match(r"(\d+)\s+\([^)]+\)\s+aligned concordantly exactly 1 time", stripped)
            if m:
                paired_aligned_one = int(m.group(1))
                continue
            m = re.match(r"(\d+)\s+\([^)]+\)\s+aligned concordantly >1 times", stripped)
            if m:
                paired_aligned_many = int(m.group(1))
                continue
            m = re.match(r"(\d+)\s+\([^)]+\)\s+aligned exactly 1 time", stripped)
            if m:
                single_aligned_one += int(m.group(1))
                continue
            m = re.match(r"(\d+)\s+\([^)]+\)\s+aligned >1 times", stripped)
            if m:
                single_aligned_many += int(m.group(1))
                continue
            m = re.match(r"(\d+)\s+\([^)]+\)\s+aligned 0 times", stripped)
            if m:
                single_unaligned += int(m.group(1))
                continue
            m = re.match(r"(\d+)\s+reads; of these:", stripped)
            if m:
                single_total = int(m.group(1))
                continue
    if paired_total:
        aligned = paired_aligned_one + paired_aligned_many
        if not paired_unaligned:
            paired_unaligned = max(paired_total - aligned, 0)
        return {"total": paired_total, "aligned": aligned, "unaligned": paired_unaligned}
    return {
        "total": single_total,
        "aligned": single_aligned_one + single_aligned_many,
        "unaligned": single_unaligned,
    }


def split_features(value: str) -> list[str]:
    if not value:
        return []
    value = str(value).strip()
    if value in {"", "-", "NA", "None", "nan"}:
        return []
    return [part.strip() for part in re.split(r"[,\s;|]+", value) if part.strip() and part.strip() != "-"]


def normalize_ko(term: str) -> str:
    term = term.strip()
    if ":" in term:
        term = term.rsplit(":", 1)[-1]
    return term


def read_feature_matrix(path: Path) -> tuple[str, list[str], dict[str, list[float]]]:
    with path.open(newline="", errors="replace") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        first = header[0]
        samples = header[1:]
        data: dict[str, list[float]] = {}
        for row in reader:
            if not row:
                continue
            values = []
            for raw in row[1 : len(samples) + 1]:
                try:
                    values.append(float(raw))
                except ValueError:
                    values.append(0.0)
            while len(values) < len(samples):
                values.append(0.0)
            data[row[0]] = values
    return first, samples, data


def table_delimiter(path: Path) -> str:
    return "," if path.suffix.lower() == ".csv" else "\t"


def read_table_header(path: Path) -> list[str]:
    with path.open(newline="", errors="replace") as handle:
        reader = csv.reader(handle, delimiter=table_delimiter(path))
        try:
            return next(reader)
        except StopIteration:
            return []


def duplicate_values(values: list[str]) -> list[str]:
    seen: set[str] = set()
    dupes: set[str] = set()
    for value in values:
        if value in seen:
            dupes.add(value)
        seen.add(value)
    return sorted(dupes)


def read_metadata_sample_names(path: Path) -> list[str]:
    with path.open(newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter=table_delimiter(path))
        fields = list(reader.fieldnames or [])
        if "sample" not in fields and "sampleID" not in fields:
            raise RuntimeError(f"Metadata must contain sample or sampleID column: {path}")
        samples = []
        for row in reader:
            sample = (row.get("sample") or row.get("sampleID") or "").strip()
            if sample:
                samples.append(sample)
    return samples


def r_string(value: Any) -> str:
    return json.dumps(str(value), ensure_ascii=True)


def parse_positive_int_list(value: Any, default: list[int]) -> list[int]:
    if value is None or value == "":
        return default
    if isinstance(value, str):
        raw_items = re.split(r"[,;\s]+", value.strip())
    elif isinstance(value, list):
        raw_items = value
    else:
        raw_items = [value]
    sizes: list[int] = []
    for raw in raw_items:
        if raw in {"", None}:
            continue
        size = int(raw)
        if size <= 0:
            raise RuntimeError(f"Expected positive integer, got: {raw}")
        if size not in sizes:
            sizes.append(size)
    return sizes or default


def write_feature_matrix(path: Path, first_col: str, samples: list[str], data: dict[str, list[float]], integer: bool = False) -> None:
    mkdir(path.parent)
    rows = sorted(data.items(), key=lambda item: sum(item[1]), reverse=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([first_col] + samples)
        for feature, values in rows:
            if integer:
                out_values = [str(int(round(v))) for v in values]
            else:
                out_values = [f"{v:.6g}" for v in values]
            writer.writerow([feature] + out_values)


def load_eggnog_annotations(path: Path) -> dict[str, dict[str, list[str]]]:
    header: list[str] | None = None
    with path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("#query"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                break
    if not header:
        raise RuntimeError(f"Cannot find #query header in eggNOG annotations: {path}")

    annotations: dict[str, dict[str, list[str]]] = {}
    with path.open(errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(header):
                fields += [""] * (len(header) - len(fields))
            row = dict(zip(header, fields))
            query = row.get("query", "").strip()
            if not query:
                continue
            annotations[query] = {
                "KO": [normalize_ko(term) for term in split_features(row.get("KEGG_ko", ""))],
                "EC": split_features(row.get("EC", "")),
                "PFAM": split_features(row.get("PFAMs", "")),
                "KEGG_Pathway": split_features(row.get("KEGG_Pathway", "")),
            }
    return annotations


def aggregate_by_mapping(
    matrix_path: Path,
    annotations: dict[str, dict[str, list[str]]],
    annotation_key: str,
) -> tuple[list[str], dict[str, list[float]]]:
    _first, samples, matrix = read_feature_matrix(matrix_path)
    agg: dict[str, list[float]] = defaultdict(lambda: [0.0] * len(samples))
    for gene, values in matrix.items():
        for feature in annotations.get(gene, {}).get(annotation_key, []):
            if not feature:
                continue
            current = agg[feature]
            for idx, value in enumerate(values):
                current[idx] += value
    return samples, dict(agg)


def rank_values(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda idx: values[idx])
    ranks = [0.0] * len(values)
    pos = 0
    while pos < len(order):
        end = pos + 1
        while end < len(order) and values[order[end]] == values[order[pos]]:
            end += 1
        avg = (pos + 1 + end) / 2.0
        for idx in order[pos:end]:
            ranks[idx] = avg
        pos = end
    return ranks


def pearson(x: list[float], y: list[float]) -> float:
    if len(x) != len(y) or len(x) < 2:
        return 0.0
    mx = sum(x) / len(x)
    my = sum(y) / len(y)
    num = sum((a - mx) * (b - my) for a, b in zip(x, y))
    denx = sum((a - mx) ** 2 for a in x)
    deny = sum((b - my) ** 2 for b in y)
    if denx <= 0 or deny <= 0:
        return 0.0
    return num / ((denx * deny) ** 0.5)


def spearman(x: list[float], y: list[float]) -> float:
    return pearson(rank_values(x), rank_values(y))


class Pipeline:
    def __init__(self, config: dict[str, Any]):
        self.config = config
        self.outdir = Path(config["output_dir"]).expanduser().resolve()
        self.envs = {**DEFAULT_ENVS, **{
            str(key): str(value)
            for key, value in config.get("envs", {}).items()
            if str(value).strip()
        }}
        self.dbs = {**DEFAULT_DBS, **{
            str(key): str(value)
            for key, value in config.get("databases", {}).items()
            if str(value).strip()
        }}
        mkdir(self.outdir)
        self.log_path = self.outdir / "pipeline.log"
        self.state_dir = self.outdir / ".pipeline_state"
        mkdir(self.state_dir)
        self.steps: set[str] = set(config.get("steps", []))
        self.samples: list[Sample] = []

        self.trim_dir = self.outdir / "01_fastp_trim_QC"
        self.rrna_dir = self.outdir / "02_rRNA_removed"
        self.assembly_dir = self.outdir / "03_assembly_megahit"
        self.contigs_dir = self.outdir / "03b_contigs_all"
        self.quast_dir = self.outdir / "03c_quast_metaquast"
        self.transdecoder_dir = self.outdir / "04_transdecoder"
        self.cds_dir = self.outdir / "04a_all_cds"
        self.pep_dir = self.outdir / "04b_all_pep"
        self.salmon_dir = self.outdir / "06_quant_salmon"
        self.salmon_matrix_dir = self.outdir / "07_salmon_matrices"
        self.eggnog_dir = self.outdir / "08_eggnog"
        self.aggregate_tpm_dir = self.outdir / "09a_aggregate_tpm_by_eggnog"
        self.aggregate_counts_dir = self.outdir / "09b_aggregate_counts_by_eggnog"
        self.kegg_dir = self.outdir / "11_KO_to_KEGG_Pathway"
        self.comparative_dir = self.outdir / "10_comparative_analysis"
        self.rnaseq_dir = self.outdir / "13_RNAseq_one_group"
        self.amr_ko_dir = self.outdir / "14_AMR_KO"
        self.card_dir = self.outdir / "15_CARD"
        self.card_summary_dir = self.outdir / "15a_CARD_analysis"
        self.integration_dir = self.outdir / "16_integracja_CARD_KO"

    def log(self, message: str) -> None:
        stamp = time.strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{stamp}] {message}"
        print(line, flush=True)
        with self.log_path.open("a") as handle:
            handle.write(line + "\n")

    def run_cmd(
        self,
        cmd: list[str],
        label: str,
        cwd: Path | None = None,
        allow_fail: bool = False,
        input_text: str | None = None,
    ) -> int:
        self.log(f"START {label}")
        self.log("CMD " + " ".join(q(part) for part in cmd))
        if input_text:
            self.log(f"STDIN {len(input_text.splitlines())} lines")
        with self.log_path.open("a") as log_handle:
            proc = subprocess.run(
                cmd,
                cwd=str(cwd) if cwd else None,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
                input=input_text,
            )
        if proc.returncode != 0:
            self.log(f"FAILED {label} exit={proc.returncode}")
            if not allow_fail:
                raise RuntimeError(f"{label} failed, see {self.log_path}")
        else:
            self.log(f"DONE {label}")
        return proc.returncode

    def run_bash(self, env_key: str, script: str, label: str, cwd: Path | None = None, allow_fail: bool = False) -> int:
        env = self.envs.get(env_key, "").strip()
        if not env:
            raise RuntimeError(
                f"Conda environment '{env_key}' is not configured. "
                f"Set envs.{env_key} in the JSON config or METAREX_ENV_{env_key.upper()}."
            )
        cmd = ["conda", "run", "-p", env, "bash", "-lc", f"set -euo pipefail\n{script}"]
        return self.run_cmd(cmd, label, cwd=cwd, allow_fail=allow_fail)

    def run_interactive_script(
        self,
        env_key: str,
        cmd_tail: list[str],
        input_lines: list[str],
        label: str,
        allow_fail: bool = False,
    ) -> int:
        env = self.envs.get(env_key, "").strip()
        if not env:
            raise RuntimeError(
                f"Conda environment '{env_key}' is not configured. "
                f"Set envs.{env_key} in the JSON config or METAREX_ENV_{env_key.upper()}."
            )
        cmd = ["conda", "run", "-p", env, *cmd_tail]
        return self.run_cmd(cmd, label, allow_fail=allow_fail, input_text="\n".join(input_lines) + "\n")

    def rscript_bin(self) -> str:
        configured = str(self.config.get("rscript_bin") or "").strip()
        if configured:
            return configured
        if Path("/usr/bin/Rscript").exists():
            return "/usr/bin/Rscript"
        return "Rscript"

    def materialize_r_script(
        self,
        script: Path,
        assignments: dict[str, Any],
        name: str,
        heatmap_sizes: list[int] | None = None,
    ) -> Path:
        text = script.read_text(encoding="utf-8")
        missing: list[str] = []
        for variable, value in assignments.items():
            pattern = re.compile(rf"^{re.escape(variable)}\s*<-\s*readline\([^\n]*\)", re.MULTILINE)
            text, replacements = pattern.subn(f"{variable} <- {r_string(value)}", text, count=1)
            if replacements == 0:
                direct_pattern = re.compile(rf"^{re.escape(variable)}\s*<-\s*[^\n]*", re.MULTILINE)
                text, replacements = direct_pattern.subn(f"{variable} <- {r_string(value)}", text, count=1)
            if replacements == 0:
                missing.append(variable)
        if missing:
            raise RuntimeError(f"Cannot make {script} non-interactive; missing readline assignments: {', '.join(missing)}")
        text = re.sub(r"(?m)^\s*\\\s*$\n?", "", text)
        text = re.sub(
            r"estimateSizeFactors\(\s*dds\s*\)",
            'estimateSizeFactors(dds, type = "poscounts")',
            text,
        )
        text = re.sub(r"DESeq\(\s*dds\s*\)", 'DESeq(dds, sfType = "poscounts")', text)
        if heatmap_sizes is not None:
            text = text.replace(
                "pdf(file.path(outdir, \"figures\", paste0(filename_base, \".pdf\")), width = 7, height = 9)",
                "pdf(file.path(outdir, \"figures\", paste0(filename_base, \".pdf\")), width = 8.5, height = 9)",
            )
            text = text.replace("show_rownames = FALSE", "show_rownames = TRUE")
            text = text.replace("cluster_rows = TRUE", "cluster_rows = FALSE", 1)
            text = text.replace(
                "show_rownames = TRUE,\n    cluster_rows = FALSE,",
                "show_rownames = TRUE,\n    fontsize_row = ifelse(ntop <= 50, 7, 5),\n    cluster_rows = FALSE,",
            )
            replacement_calls = [f'make_heatmap({size}, "Heatmap_vst_top{size}")' for size in heatmap_sizes]
            updated_lines = []
            replaced = 0
            for line in text.splitlines():
                if re.match(r'\s*make_heatmap\(\s*\d+\s*,\s*"Heatmap_vst_top\d+"\s*\)\s*$', line):
                    if replaced == 0:
                        updated_lines.extend(replacement_calls)
                    replaced += 1
                else:
                    updated_lines.append(line)
            if replaced == 0:
                raise RuntimeError(f"Cannot configure heatmap sizes in {script}; make_heatmap calls were not found")
            text = "\n".join(updated_lines) + "\n"
        out = self.state_dir / "r_scripts"
        mkdir(out)
        safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_") or script.stem
        materialized = out / f"{safe_name}.R"
        materialized.write_text(text, encoding="utf-8")
        return materialized

    def run_r_script(
        self,
        script: Path,
        label: str,
        allow_fail: bool = False,
    ) -> int:
        return self.run_cmd(
            [self.rscript_bin(), str(script)],
            label,
            allow_fail=allow_fail,
        )

    def planned_step_names(self) -> set[str]:
        start_from = self.config.get("start_from", "beginning")
        started = start_from in {"", "beginning"}
        planned: set[str] = set()
        for step_name, _label in PIPELINE_STEPS:
            if step_name == start_from:
                started = True
            if started and step_name in self.steps:
                planned.add(step_name)
        return planned

    def run(self) -> None:
        self.log("RNA-seq AMR pipeline started")
        self.log("Source script catalogue: " + json.dumps(SOURCE_SCRIPTS, sort_keys=True))
        input_dir = Path(self.config["input_dir"]).expanduser().resolve()
        recursive = bool(self.config.get("recursive", True))
        planned_steps = self.planned_step_names()

        existing_samples = load_sample_table(self.outdir / "00_samples.tsv")
        materialize_steps = {
            "estimate_rrna",
            "trim_qc",
            "remove_rrna",
            "megahit",
            "salmon",
            "card_rgi",
        }
        self.samples = discover_samples(input_dir, recursive=recursive)
        if existing_samples and not planned_steps & materialize_steps:
            self.samples = existing_samples
            self.log("Using existing 00_samples.tsv; selected steps do not need raw-read materialization")
        elif not self.samples:
            if existing_samples:
                self.samples = existing_samples
                self.log("Using existing 00_samples.tsv because no new FASTQ pairs were discovered")
            else:
                raise RuntimeError(f"No paired FASTQ samples found in {input_dir}")

        if self.config.get("run_profile", "full") == "test" and planned_steps & materialize_steps:
            self.apply_test_run_if_requested()

        sample_table(self.samples, self.outdir)
        self.write_metadata()
        self.log(f"Discovered {len(self.samples)} paired-end samples")

        ordered_steps: list[tuple[str, Any]] = [
            ("estimate_rrna", self.estimate_rrna_content),
            ("trim_qc", self.trim_qc),
            ("remove_rrna", self.remove_rrna),
            ("rrna_stats", self.rrna_stats),
            ("megahit", self.megahit),
            ("collect_contigs", self.collect_contigs),
            ("quast", self.quast),
            ("transdecoder", self.transdecoder),
            ("collect_cds", self.collect_cds),
            ("cds_stats", self.cds_stats),
            ("salmon", self.salmon),
            ("salmon_matrix", self.salmon_matrices),
            ("eggnog", self.eggnog),
            ("aggregate", self.aggregate_annotations),
            ("kegg_pathways", self.kegg_pathways),
            ("comparative", self.comparative_analysis),
            ("rnaseq_overview", self.rnaseq_overview),
            ("amr_ko", self.amr_ko),
            ("card_rgi", self.card_rgi),
            ("card_summary", self.card_summary),
            ("card_ko_integration", self.card_ko_integration),
            ("multiqc", self.final_multiqc),
            ("summary", self.report),
        ]
        start_from = self.config.get("start_from", "beginning")
        started = start_from in {"", "beginning"}
        for step_name, func in ordered_steps:
            if step_name == start_from:
                started = True
            if not started:
                self.log(f"SKIP {step_name}: before start_from={start_from}")
                continue
            if step_name not in self.steps:
                continue
            self.run_step(step_name, func)

        if "summary" not in self.steps:
            self.report()
        self.log("RNA-seq AMR pipeline finished")

    def marker(self, step_name: str, suffix: str) -> Path:
        return self.state_dir / f"{step_name}.{suffix}"

    def write_step_state(self, step_name: str, status: str, error: str = "") -> None:
        payload = {
            "step": step_name,
            "status": status,
            "time": time.strftime("%Y-%m-%d %H:%M:%S"),
            "error": error,
        }
        (self.state_dir / f"{step_name}.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    def run_step(self, step_name: str, func: Any) -> None:
        done = self.marker(step_name, "done")
        failed = self.marker(step_name, "failed")
        resume = bool(self.config.get("resume", True))
        force = bool(self.config.get("force", False))
        if resume and not force and done.exists():
            if self.step_outputs_present(step_name):
                self.log(f"SKIP {step_name}: checkpoint exists ({done})")
                return
            self.log(f"RERUN {step_name}: checkpoint exists but expected outputs are missing")

        failed.unlink(missing_ok=True)
        self.marker(step_name, "running").write_text(time.strftime("%Y-%m-%d %H:%M:%S") + "\n", encoding="utf-8")
        self.write_step_state(step_name, "running")
        self.log(f"CHECKPOINT START {step_name}")
        try:
            func()
        except Exception as exc:
            self.marker(step_name, "running").unlink(missing_ok=True)
            failed.write_text(f"{time.strftime('%Y-%m-%d %H:%M:%S')}\n{exc}\n", encoding="utf-8")
            self.write_step_state(step_name, "failed", str(exc))
            self.log(f"CHECKPOINT FAILED {step_name}: {exc}")
            raise
        self.marker(step_name, "running").unlink(missing_ok=True)
        done.write_text(time.strftime("%Y-%m-%d %H:%M:%S") + "\n", encoding="utf-8")
        self.write_step_state(step_name, "done")
        self.log(f"CHECKPOINT DONE {step_name}")

    def step_outputs_present(self, step_name: str) -> bool:
        if step_name == "estimate_rrna":
            return (self.outdir / "00_rRNA_content" / "rrna_content_summary.tsv").exists()
        if step_name == "trim_qc":
            clean = self.trim_dir / "fastq_clean"
            return all((clean / f"{s.sample}_R1_clean.fastq.gz").exists() and (clean / f"{s.sample}_R2_clean.fastq.gz").exists() for s in self.samples)
        if step_name == "remove_rrna":
            return all((self.rrna_dir / f"{s.sample}.rRNAfree_R1.fastq.gz").exists() and (self.rrna_dir / f"{s.sample}.rRNAfree_R2.fastq.gz").exists() for s in self.samples)
        if step_name == "rrna_stats":
            return (self.outdir / "02a_rRNA_stats" / "rrna_removal_stats.tsv").exists()
        if step_name == "megahit":
            return all((self.assembly_dir / s.sample / "final.contigs.fa").exists() for s in self.samples)
        if step_name == "collect_contigs":
            return bool(list(self.contigs_dir.glob("*.fa")))
        if step_name == "quast":
            return (self.quast_dir / "report.tsv").exists()
        if step_name == "transdecoder":
            return bool(list(self.transdecoder_dir.glob("*/*.transdecoder.cds")))
        if step_name == "collect_cds":
            return bool(list(self.cds_dir.glob("*.cds.fa")))
        if step_name == "cds_stats":
            return (self.outdir / "05_cds_stats" / "cds_comparison.tsv").exists()
        if step_name == "salmon":
            return bool(list(self.salmon_dir.glob("*/quant.sf")))
        if step_name == "salmon_matrix":
            return (self.salmon_matrix_dir / "salmon_TPM_all.tsv").exists()
        if step_name == "eggnog":
            return (self.eggnog_dir / "eggnog_cds.emapper.annotations").exists()
        if step_name == "aggregate":
            return (self.aggregate_counts_dir / "09b_COUNTS_KO_matrix.tsv").exists()
        if step_name == "kegg_pathways":
            return (self.kegg_dir / "11_COUNTS_KEGG_Pathway_matrix.tsv").exists()
        if step_name == "comparative":
            return all(
                path.exists()
                for path in [
                    self.comparative_dir / "10a_DEG_KO_ALDEx2" / "10_ALDEx2_KO_all_comparisons.tsv",
                    self.comparative_dir / "10b_DEG_KO_DESeq2" / "DESeq2_KO_all_comparisons.tsv",
                    self.comparative_dir / "10a_DEG_KEGG_Pathway_ALDEx2" / "ALDEx2_pairwise_all_comparisons.tsv",
                    self.comparative_dir / "10b_DEG_KEGG_Pathway_DESeq2" / "DESeq2_KEGG_Pathway_all_comparisons.tsv",
                ]
            )
        if step_name == "rnaseq_overview":
            return (self.rnaseq_dir / "tables" / "PCA_coordinates.tsv").exists()
        if step_name == "amr_ko":
            return (self.amr_ko_dir / "tables" / "AMR_activity_index.tsv").exists()
        if step_name == "card_rgi":
            return bool(list(self.card_dir.glob("*/card_amr_*.txt")))
        if step_name == "card_summary":
            return (self.card_summary_dir / "tables" / "CARD_all_filtered.tsv").exists()
        if step_name == "card_ko_integration":
            return (self.integration_dir / "tables" / "CARD_KO_integrated_metrics.tsv").exists()
        if step_name == "multiqc":
            return (self.outdir / "99_multiqc_all" / "multiqc_all.html").exists()
        if step_name == "summary":
            return (self.outdir / "index.html").exists()
        return False

    def apply_test_run_if_requested(self) -> None:
        sample_count = max(1, int(self.config.get("test_sample_count", 1)))
        read_count = max(1, int(self.config.get("test_read_count", 10000)))
        seed = int(self.config.get("test_seed", 7))
        selected = self.samples[:sample_count]
        out = self.outdir / "00_test_run_downsampled"
        mkdir(out)
        rows = [["sample", "input_r1", "input_r2", "output_r1", "output_r2", "reads_written"]]
        test_samples: list[Sample] = []
        self.log(f"TEST RUN enabled: samples={len(selected)}/{len(self.samples)}, target_reads={read_count}, seed={seed}")
        for idx, sample in enumerate(selected, 1):
            if not sample.r1 or not sample.r2:
                continue
            out_r1 = out / f"{sample.sample}_R1.test.fastq.gz"
            out_r2 = out / f"{sample.sample}_R2.test.fastq.gz"
            if out_r1.exists() and out_r2.exists() and not bool(self.config.get("force", False)):
                written = "existing"
            else:
                written = str(write_random_paired_fastq_subset(sample.r1, sample.r2, out_r1, out_r2, read_count, seed + idx))
            test_samples.append(Sample(sample=sample.sample, r1=out_r1, r2=out_r2))
            rows.append([sample.sample, sample.r1, sample.r2, out_r1, out_r2, written])
        with (out / "manifest.tsv").open("w", newline="") as handle:
            csv.writer(handle, delimiter="\t").writerows(rows)
        self.samples = test_samples

    def write_metadata(self) -> Path:
        meta_out = self.outdir / "00_metadata.tsv"
        metadata_file = self.config.get("metadata_file", "").strip()
        sample_names = [sample.sample for sample in self.samples]
        if metadata_file:
            src = Path(metadata_file).expanduser()
            if not src.exists():
                raise RuntimeError(f"Metadata file does not exist: {src}")
            delimiter = table_delimiter(src)
            with src.open(newline="", errors="replace") as handle:
                reader = csv.DictReader(handle, delimiter=delimiter)
                rows = list(reader)
                fields = list(reader.fieldnames or [])
            if "sample" not in fields and "sampleID" not in fields:
                raise RuntimeError("Metadata file must contain a sample or sampleID column")
            if "sample" not in fields and "sampleID" in fields:
                fields.append("sample")
                for row in rows:
                    row["sample"] = row.get("sampleID", "")
            if "sampleID" not in fields and "sample" in fields:
                fields.append("sampleID")
                for row in rows:
                    row["sampleID"] = row.get("sample", "")
            if "group" not in fields:
                fields.append("group")
                for row in rows:
                    row["group"] = "all"
            sample_set = set(sample_names)
            rows_by_sample: dict[str, dict[str, str]] = {}
            ignored = 0
            for row in rows:
                sample = (row.get("sample") or row.get("sampleID") or "").strip()
                if not sample or sample not in sample_set:
                    ignored += 1
                    continue
                if sample in rows_by_sample:
                    continue
                row["sample"] = sample
                row["sampleID"] = row.get("sampleID", "") or sample
                row["group"] = row.get("group", "") or "all"
                rows_by_sample[sample] = row
            output_rows = []
            missing = []
            for sample in sample_names:
                if sample in rows_by_sample:
                    output_rows.append(rows_by_sample[sample])
                    continue
                missing.append(sample)
                row = {field: "" for field in fields}
                row["sample"] = sample
                row["sampleID"] = sample
                row["group"] = "all"
                output_rows.append(row)
            with meta_out.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
                writer.writeheader()
                writer.writerows(output_rows)
            if ignored:
                self.log(f"Metadata: ignored {ignored} rows not present in current sample set")
            if missing:
                self.log("Metadata: added default rows for missing samples: " + ", ".join(missing))
            return meta_out

        with meta_out.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "sampleID", "group"])
            for sample in sample_names:
                writer.writerow([sample, sample, "all"])
        return meta_out

    def estimate_rrna_content(self) -> None:
        idx = self.config.get("rrna_index") or self.dbs["rrna_bowtie2"]
        if not idx:
            raise RuntimeError("rRNA Bowtie2 index is not configured; set rrna_index or databases.rrna_bowtie2")
        if not bowtie2_index_exists(idx):
            raise RuntimeError(f"Bowtie2 rRNA index not found for prefix: {idx}")
        out = self.outdir / "00_rRNA_content"
        logdir = out / "logs"
        mkdir(logdir)
        threads = int(self.config.get("threads", 8))
        summary = out / "rrna_content_summary.tsv"
        rows = [["Sample", "Total_Reads", "rRNA_Aligned", "mRNA_Unaligned", "Percent_rRNA", "Percent_mRNA"]]
        for sample in self.samples:
            if not sample.r1 or not sample.r2:
                continue
            log = logdir / f"bowtie2_{sample.sample}.log"
            script = (
                f"bowtie2 -x {q(idx)} -1 {q(sample.r1)} -2 {q(sample.r2)} "
                f"--end-to-end --very-sensitive -L 31 -N 0 --score-min L,0,-0.6 "
                f"-p {threads} -S /dev/null 2> {q(log)}"
            )
            self.run_bash("qc", script, f"Bowtie2 rRNA estimate {sample.sample}")
            stats = parse_bowtie2_log(log)
            total = stats["total"]
            aligned = stats["aligned"]
            unaligned = stats["unaligned"]
            pct_rrna = (aligned / total * 100.0) if total else 0.0
            pct_mrna = (unaligned / total * 100.0) if total else 0.0
            rows.append([sample.sample, total, aligned, unaligned, f"{pct_rrna:.2f}", f"{pct_mrna:.2f}"])
        with summary.open("w", newline="") as handle:
            csv.writer(handle, delimiter="\t").writerows(rows)

    def fastp_options(self) -> str:
        strategy = self.config.get("fastp_strategy", "full")
        phred = int(self.config.get("fastp_phred", 20))
        min_len = int(self.config.get("fastp_min_len", 50))
        presets = {
            "adapters": "--detect_adapter_for_pe",
            "ends": "--detect_adapter_for_pe --cut_front --cut_tail",
            "quality": f"--detect_adapter_for_pe --qualified_quality_phred {phred} --length_required {min_len}",
            "full": f"--detect_adapter_for_pe --cut_front --cut_tail --qualified_quality_phred {phred} --length_required {min_len}",
            "custom": self.config.get("fastp_extra", "--detect_adapter_for_pe"),
        }
        return presets.get(strategy, presets["full"])

    def trim_qc(self) -> None:
        clean = self.trim_dir / "fastq_clean"
        reports = self.trim_dir / "reports"
        mkdir(clean)
        mkdir(reports)
        threads = int(self.config.get("threads", 8))
        opts = self.fastp_options()
        for sample in self.samples:
            if not sample.r1 or not sample.r2:
                continue
            r1_out = clean / f"{sample.sample}_R1_clean.fastq.gz"
            r2_out = clean / f"{sample.sample}_R2_clean.fastq.gz"
            script = (
                f"fastp -i {q(sample.r1)} -I {q(sample.r2)} "
                f"-o {q(r1_out)} -O {q(r2_out)} "
                f"-h {q(reports / (sample.sample + '_fastp.html'))} "
                f"-j {q(reports / (sample.sample + '_fastp.json'))} "
                f"-w {threads} {opts}"
            )
            self.run_bash("qc", script, f"fastp {sample.sample}")
        fastp_inputs = " ".join(q(path) for path in reports.glob("*_fastp.json"))
        if fastp_inputs:
            self.run_bash("qc", f"multiqc --force --module fastp {fastp_inputs} -o {q(self.trim_dir / 'multiqc')} -n multiqc_fastp.html", "MultiQC fastp", allow_fail=True)

    def clean_reads_for_sample(self, sample: Sample) -> tuple[Path, Path]:
        clean = self.trim_dir / "fastq_clean"
        r1 = clean / f"{sample.sample}_R1_clean.fastq.gz"
        r2 = clean / f"{sample.sample}_R2_clean.fastq.gz"
        if r1.exists() and r2.exists():
            return r1, r2
        if sample.r1 and sample.r2:
            return sample.r1, sample.r2
        raise RuntimeError(f"No paired reads available for sample {sample.sample}")

    def rrna_free_reads_for_sample(self, sample: Sample) -> tuple[Path, Path] | None:
        r1 = self.rrna_dir / f"{sample.sample}.rRNAfree_R1.fastq.gz"
        r2 = self.rrna_dir / f"{sample.sample}.rRNAfree_R2.fastq.gz"
        if r1.exists() and r2.exists():
            return r1, r2
        return None

    def assembly_reads_for_sample(self, sample: Sample) -> tuple[Path, Path]:
        rrna = self.rrna_free_reads_for_sample(sample)
        if rrna:
            return rrna
        return self.clean_reads_for_sample(sample)

    def remove_rrna(self) -> None:
        idx = self.config.get("rrna_index") or self.dbs["rrna_bowtie2"]
        if not idx:
            raise RuntimeError("rRNA Bowtie2 index is not configured; set rrna_index or databases.rrna_bowtie2")
        if not bowtie2_index_exists(idx):
            raise RuntimeError(f"Bowtie2 rRNA index not found for prefix: {idx}")
        mkdir(self.rrna_dir)
        logdir = self.rrna_dir / "bowtie2_logs"
        mkdir(logdir)
        threads = int(self.config.get("threads", 8))
        for sample in self.samples:
            r1, r2 = self.clean_reads_for_sample(sample)
            out1 = self.rrna_dir / f"{sample.sample}.rRNAfree_R1.fastq.gz"
            out2 = self.rrna_dir / f"{sample.sample}.rRNAfree_R2.fastq.gz"
            tmp_pattern = self.rrna_dir / f"{sample.sample}.rRNAfree_%.fastq.gz"
            tmp1 = self.rrna_dir / f"{sample.sample}.rRNAfree_1.fastq.gz"
            tmp2 = self.rrna_dir / f"{sample.sample}.rRNAfree_2.fastq.gz"
            log = logdir / f"{sample.sample}.log"
            script = (
                f"rm -f {q(tmp1)} {q(tmp2)} {q(out1)} {q(out2)}\n"
                f"bowtie2 -x {q(idx)} -1 {q(r1)} -2 {q(r2)} "
                f"--end-to-end --very-sensitive -L 31 -N 0 --score-min L,0,-0.6 "
                f"--threads {threads} --un-conc-gz {q(tmp_pattern)} -S /dev/null 2> {q(log)}\n"
                f"test -f {q(tmp1)} && test -f {q(tmp2)}\n"
                f"mv -f {q(tmp1)} {q(out1)}\n"
                f"mv -f {q(tmp2)} {q(out2)}"
            )
            self.run_bash("qc", script, f"Bowtie2 remove rRNA {sample.sample}")

    def rrna_stats(self) -> None:
        out = self.outdir / "02a_rRNA_stats"
        mkdir(out)
        rows = [["Sample", "Reads_input", "Reads_clean", "Reads_rRNAfree", "Percent_remaining"]]
        missing: list[str] = []
        for sample in self.samples:
            if not sample.r1 or not sample.r2:
                continue
            try:
                clean1, clean2 = self.clean_reads_for_sample(sample)
                rrna = self.rrna_free_reads_for_sample(sample)
                if not rrna:
                    raise RuntimeError("missing rRNA-free FASTQ")
                input_reads = count_fastq_reads(sample.r1) + count_fastq_reads(sample.r2)
                clean_reads = count_fastq_reads(clean1) + count_fastq_reads(clean2)
                rrna_reads = count_fastq_reads(rrna[0]) + count_fastq_reads(rrna[1])
                pct = (rrna_reads / clean_reads * 100.0) if clean_reads else 0.0
                rows.append([sample.sample, input_reads, clean_reads, rrna_reads, f"{pct:.2f}"])
            except Exception as exc:
                missing.append(f"{sample.sample}\t{exc}")
        with (out / "rrna_removal_stats.tsv").open("w", newline="") as handle:
            csv.writer(handle, delimiter="\t").writerows(rows)
        if missing:
            (out / "missing_files.log").write_text("\n".join(missing) + "\n", encoding="utf-8")

    def megahit(self) -> None:
        mkdir(self.assembly_dir)
        threads = int(self.config.get("threads", 16))
        klist = self.config.get("megahit_k_list", "21,41,61,81,101,121")
        min_len = int(self.config.get("megahit_min_contig_len", 300))
        failures = 0
        for sample in self.samples:
            r1, r2 = self.assembly_reads_for_sample(sample)
            sample_out = self.assembly_dir / sample.sample
            if (sample_out / "final.contigs.fa").exists() and not bool(self.config.get("force", False)):
                self.log(f"SKIP MEGAHIT {sample.sample}: final.contigs.fa exists")
                continue
            script = f"megahit -1 {q(r1)} -2 {q(r2)} -o {q(sample_out)} -t {threads} --k-list {q(klist)} --min-contig-len {min_len}"
            rc = self.run_bash("qc", script, f"MEGAHIT {sample.sample}", allow_fail=True)
            if rc != 0 or not (sample_out / "final.contigs.fa").exists():
                failures += 1
        if failures:
            raise RuntimeError(f"MEGAHIT failed for {failures}/{len(self.samples)} samples")

    def collect_contigs(self) -> None:
        mkdir(self.contigs_dir)
        copied = 0
        for sample in self.samples:
            src = self.assembly_dir / sample.sample / "final.contigs.fa"
            if not src.exists():
                continue
            dst = self.contigs_dir / f"{sample.sample}.contigs.fa"
            if not dst.exists() or dst.stat().st_mtime < src.stat().st_mtime:
                shutil.copy2(src, dst)
            copied += 1
        self.log(f"Collected {copied} contig FASTA files in {self.contigs_dir}")

    def quast(self) -> None:
        files = sorted(self.contigs_dir.glob("*.fa"))
        if not files:
            self.collect_contigs()
            files = sorted(self.contigs_dir.glob("*.fa"))
        if not files:
            raise RuntimeError("No contigs available for QUAST")
        mkdir(self.quast_dir)
        threads = int(self.config.get("threads", 16))
        quast_mode = self.config.get("quast_mode", "quast")
        quast_bin = "metaquast.py" if quast_mode == "metaquast" else "quast.py"
        no_reference_args = "--max-ref-number 0 " if quast_mode == "metaquast" else ""
        labels = ",".join(path.stem.replace(".contigs", "") for path in files)
        script = f"{quast_bin} -t {threads} --min-contig 200 {no_reference_args}-o {q(self.quast_dir)} --labels {q(labels)} " + " ".join(q(path) for path in files)
        rc = self.run_bash("quast", script, "QUAST/MetaQUAST de novo report", allow_fail=True)
        if rc != 0 or not (self.quast_dir / "report.tsv").exists():
            raise RuntimeError(f"QUAST/MetaQUAST failed or report.tsv missing in {self.quast_dir}")
        self.run_bash("qc", f"multiqc --force --module quast --ignore '*multiqc*' {q(self.quast_dir)} -o {q(self.quast_dir)} -n multiqc_quast.html", "MultiQC QUAST", allow_fail=True)

    def transdecoder(self) -> None:
        files = sorted(self.contigs_dir.glob("*.fa"))
        if not files:
            self.collect_contigs()
            files = sorted(self.contigs_dir.glob("*.fa"))
        if not files:
            raise RuntimeError("No contigs available for TransDecoder")
        for fasta in files:
            sample = fasta.name
            sample = re.sub(r"\.contigs\.fa$", "", sample)
            sample = re.sub(r"\.fa$", "", sample)
            workdir = self.transdecoder_dir / sample
            mkdir(workdir)
            transcript = workdir / f"{sample}.transcripts.fa"
            if not transcript.exists() or transcript.stat().st_mtime < fasta.stat().st_mtime:
                shutil.copy2(fasta, transcript)
            script = (
                f"TransDecoder.LongOrfs -t {q(transcript.name)}\n"
                f"TransDecoder.Predict -t {q(transcript.name)} --no_refine_starts"
            )
            self.run_bash("meta", script, f"TransDecoder {sample}", cwd=workdir)

    def collect_cds(self) -> None:
        mkdir(self.cds_dir)
        mkdir(self.pep_dir)
        cds_files = sorted(self.transdecoder_dir.glob("*/*.transdecoder.cds"))
        pep_files = sorted(self.transdecoder_dir.glob("*/*.transdecoder.pep"))
        if not cds_files:
            raise RuntimeError(f"No *.transdecoder.cds files found in {self.transdecoder_dir}")
        for src in cds_files:
            sample = src.parent.name
            copy_fasta_with_sample_prefix(src, self.cds_dir / f"{sample}.cds.fa", sample)
        for src in pep_files:
            sample = src.parent.name
            copy_fasta_with_sample_prefix(src, self.pep_dir / f"{sample}.transdecoder.pep", sample)
        self.log(f"Collected {len(cds_files)} CDS and {len(pep_files)} protein FASTA files")

    def cds_stats(self) -> None:
        out = self.outdir / "05_cds_stats"
        mkdir(out)
        rows = [["Sample", "CDS_count", "Min_length", "Max_length", "Mean_length", "Median_length", "Total_length"]]
        for cds in sorted(self.cds_dir.glob("*.cds.fa")):
            lengths = fasta_lengths(cds)
            sample = cds.name.replace(".cds.fa", "")
            if lengths:
                rows.append([
                    sample,
                    len(lengths),
                    min(lengths),
                    max(lengths),
                    f"{sum(lengths) / len(lengths):.2f}",
                    f"{median(lengths):.2f}",
                    sum(lengths),
                ])
            else:
                rows.append([sample, 0, 0, 0, "0.00", "0.00", 0])
        with (out / "cds_comparison.tsv").open("w", newline="") as handle:
            csv.writer(handle, delimiter="\t").writerows(rows)

    def salmon(self) -> None:
        cds_files = sorted(self.cds_dir.glob("*.cds.fa"))
        if not cds_files:
            raise RuntimeError(f"No CDS FASTA files found in {self.cds_dir}")
        threads = int(self.config.get("salmon_threads", self.config.get("threads", 8)))
        sample_map = {sample.sample: sample for sample in self.samples}
        failures = 0
        for cds in cds_files:
            sample_name = cds.name.replace(".cds.fa", "")
            sample = sample_map.get(sample_name)
            if not sample:
                self.log(f"SKIP Salmon {sample_name}: sample not in sample table")
                continue
            r1, r2 = self.assembly_reads_for_sample(sample)
            workdir = self.salmon_dir / sample_name
            mkdir(workdir)
            script = (
                f"salmon index -t {q(cds)} -i {q(workdir / 'salmon_index')}\n"
                f"salmon quant -i {q(workdir / 'salmon_index')} -l A -1 {q(r1)} -2 {q(r2)} -p {threads} -o {q(workdir)}"
            )
            rc = self.run_bash("meta", script, f"Salmon {sample_name}", allow_fail=True)
            if rc != 0 or not (workdir / "quant.sf").exists():
                failures += 1
        if failures:
            raise RuntimeError(f"Salmon failed for {failures}/{len(cds_files)} CDS files")

    def salmon_matrices(self) -> None:
        quant_files = sorted(self.salmon_dir.glob("*/quant.sf"))
        if not quant_files:
            raise RuntimeError(f"No quant.sf files found in {self.salmon_dir}")
        tpm: dict[str, dict[str, float]] = {}
        counts: dict[str, dict[str, float]] = {}
        genes: set[str] = set()
        samples: list[str] = []
        for qf in quant_files:
            sample = qf.parent.name
            samples.append(sample)
            with qf.open(newline="", errors="replace") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    gene = row.get("Name", "")
                    if not gene:
                        continue
                    genes.add(gene)
                    tpm.setdefault(gene, {})[sample] = float(row.get("TPM", 0) or 0)
                    counts.setdefault(gene, {})[sample] = float(row.get("NumReads", 0) or 0)

        mkdir(self.salmon_matrix_dir)
        tpm_matrix = {gene: [tpm.get(gene, {}).get(sample, 0.0) for sample in samples] for gene in sorted(genes)}
        counts_matrix = {gene: [counts.get(gene, {}).get(sample, 0.0) for sample in samples] for gene in sorted(genes)}
        write_feature_matrix(self.salmon_matrix_dir / "salmon_TPM_all.tsv", "Name", samples, tpm_matrix)
        write_feature_matrix(self.salmon_matrix_dir / "salmon_COUNTS_CDS_matrix.tsv", "CDS", samples, counts_matrix, integer=True)

        threshold = float(self.config.get("tpm_threshold", 1.0))
        min_samples = int(self.config.get("tpm_min_samples", 1))
        filtered = {
            gene: values
            for gene, values in tpm_matrix.items()
            if sum(1 for value in values if value >= threshold) >= min_samples
        }
        safe_threshold = str(threshold).replace(".", "_")
        write_feature_matrix(
            self.salmon_matrix_dir / f"salmon_TPM_TPM{safe_threshold}_samples{min_samples}.tsv",
            "Name",
            samples,
            filtered,
        )
        self.log(f"Salmon matrices: {len(samples)} samples, {len(genes)} CDS, {len(filtered)} passing TPM filter")

    def eggnog(self) -> None:
        pep_files = sorted(self.pep_dir.glob("*.transdecoder.pep"))
        if not pep_files:
            raise RuntimeError(f"No protein FASTA files found in {self.pep_dir}")
        mkdir(self.eggnog_dir)
        all_pep = self.eggnog_dir / "all_cds.pep"
        with all_pep.open("w") as out:
            for pep in pep_files:
                with pep.open(errors="replace") as inp:
                    shutil.copyfileobj(inp, out)
        data_dir = self.config.get("eggnog_data_dir") or self.dbs["eggnog_data"]
        if not data_dir:
            raise RuntimeError("eggNOG data directory is not configured; set eggnog_data_dir or databases.eggnog_data")
        threads = int(self.config.get("threads", 16))
        script = (
            f"export EGGNOG_DATA_DIR={q(data_dir)}\n"
            f"mkdir -p {q(self.eggnog_dir / 'tmp')}\n"
            f"emapper.py -i {q(all_pep)} --itype proteins -o eggnog_cds --cpu {threads} "
            f"--output_dir {q(self.eggnog_dir)} --temp_dir {q(self.eggnog_dir / 'tmp')} "
            f"--data_dir {q(data_dir)} --go_evidence non-electronic --override"
        )
        self.run_bash("eggnog", script, "eggNOG-mapper")

    def aggregate_annotations(self) -> None:
        anno = self.eggnog_dir / "eggnog_cds.emapper.annotations"
        if not anno.exists():
            raise RuntimeError(f"eggNOG annotations not found: {anno}")
        annotations = load_eggnog_annotations(anno)

        tpm_all = self.salmon_matrix_dir / "salmon_TPM_all.tsv"
        counts_all = self.salmon_matrix_dir / "salmon_COUNTS_CDS_matrix.tsv"
        if not tpm_all.exists() or not counts_all.exists():
            raise RuntimeError("Salmon matrices are missing; run salmon_matrix first")

        for key, label, first_col in [("KO", "KO", "KO"), ("EC", "EC", "EC"), ("PFAM", "PFAM", "PFAM")]:
            samples, tpm_agg = aggregate_by_mapping(tpm_all, annotations, key)
            write_feature_matrix(self.aggregate_tpm_dir / f"09_TPM_{label}_matrix.tsv", first_col, samples, tpm_agg)
            samples_counts, count_agg = aggregate_by_mapping(counts_all, annotations, key)
            write_feature_matrix(self.aggregate_counts_dir / f"09b_COUNTS_{label}_matrix.tsv", first_col, samples_counts, count_agg, integer=True)

        first, samples, counts_matrix = read_feature_matrix(counts_all)
        write_feature_matrix(self.aggregate_counts_dir / "09b_COUNTS_CDS_matrix.tsv", first, samples, counts_matrix, integer=True)
        self.log(f"Aggregated eggNOG annotations for {len(annotations)} CDS/proteins")

    def kegg_pathways(self) -> None:
        anno = self.eggnog_dir / "eggnog_cds.emapper.annotations"
        ko_matrix_path = self.aggregate_counts_dir / "09b_COUNTS_KO_matrix.tsv"
        if not anno.exists() or not ko_matrix_path.exists():
            raise RuntimeError("Missing eggNOG annotations or KO count matrix")
        annotations = load_eggnog_annotations(anno)
        ko_to_pathway: dict[str, set[str]] = defaultdict(set)
        for row in annotations.values():
            for ko in row.get("KO", []):
                for pathway in row.get("KEGG_Pathway", []):
                    ko_to_pathway[ko].add(pathway)
        _first, samples, ko_matrix = read_feature_matrix(ko_matrix_path)
        pathway_matrix: dict[str, list[float]] = defaultdict(lambda: [0.0] * len(samples))
        for ko, values in ko_matrix.items():
            for pathway in ko_to_pathway.get(normalize_ko(ko), set()):
                current = pathway_matrix[pathway]
                for idx, value in enumerate(values):
                    current[idx] += value
        write_feature_matrix(self.kegg_dir / "11_COUNTS_KEGG_Pathway_matrix.tsv", "KEGG_Pathway", samples, dict(pathway_matrix), integer=True)
        self.log(f"Aggregated {len(pathway_matrix)} KEGG pathways from KO counts")

    def comparative_metadata(self, group_col: str, matrix_samples: list[str]) -> tuple[Path, str]:
        meta_path = self.outdir / "00_metadata.tsv"
        if not meta_path.exists():
            raise RuntimeError(f"Metadata missing: {meta_path}")
        with meta_path.open(newline="", errors="replace") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            rows = list(reader)
            fields = list(reader.fieldnames or [])
        if group_col not in fields:
            raise RuntimeError(f"Comparative analysis requires metadata column '{group_col}'")
        sample_field = "sampleID" if "sampleID" in fields else "sample"
        if sample_field not in fields:
            raise RuntimeError("Metadata must contain sampleID or sample column for comparative analysis")
        rows_by_sample = {(row.get(sample_field) or "").strip(): row for row in rows}
        group_levels: list[str] = []
        output_rows: list[dict[str, str]] = []
        for sample in matrix_samples:
            row = rows_by_sample.get(sample)
            if not row:
                raise RuntimeError(f"Sample missing from metadata for comparative analysis: {sample}")
            group = (row.get(group_col) or "").strip()
            if not group:
                raise RuntimeError(f"Metadata column '{group_col}' is empty for sample: {sample}")
            if group not in group_levels:
                group_levels.append(group)
            out_row = {
                "sampleID": sample,
                "sample": sample,
                "group": group,
                "condition": group,
            }
            if "batch" in fields:
                out_row["batch"] = row.get("batch", "")
            output_rows.append(out_row)
        if len(group_levels) < 2:
            raise RuntimeError(f"Comparative analysis needs at least two groups in '{group_col}', got: {', '.join(group_levels) or 'none'}")
        ref_level = str(self.config.get("comparative_ref_level") or "").strip() or group_levels[0]
        if ref_level not in group_levels:
            raise RuntimeError(f"comparative_ref_level={ref_level} is not present in metadata column '{group_col}'")
        mkdir(self.comparative_dir)
        out = self.comparative_dir / f"metadata_{group_col}.tsv"
        fieldnames = ["sampleID", "sample", "group", "condition"] + (["batch"] if any("batch" in row for row in output_rows) else [])
        with out.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(output_rows)
        counts = {group: sum(1 for row in output_rows if row["group"] == group) for group in group_levels}
        self.log(f"Comparative metadata: group_col={group_col}, ref_level={ref_level}, groups={counts}")
        return out, ref_level

    def comparative_analysis(self) -> None:
        ko_matrix = self.aggregate_counts_dir / "09b_COUNTS_KO_matrix.tsv"
        pathway_matrix = self.kegg_dir / "11_COUNTS_KEGG_Pathway_matrix.tsv"
        if not ko_matrix.exists():
            raise RuntimeError(f"KO count matrix missing: {ko_matrix}")
        if not pathway_matrix.exists():
            raise RuntimeError(f"KEGG Pathway count matrix missing: {pathway_matrix}; run kegg_pathways first")
        _first, matrix_samples, _matrix = read_feature_matrix(ko_matrix)
        group_col = str(self.config.get("comparative_group_col") or "group_all").strip() or "group_all"
        meta_file, ref_level = self.comparative_metadata(group_col, matrix_samples)

        jobs = [
            (
                "compare_ko_aldex2",
                self.comparative_dir / "10a_DEG_KO_ALDEx2",
                {
                    "ko_file": ko_matrix,
                    "meta_file": meta_file,
                    "out_dir": self.comparative_dir / "10a_DEG_KO_ALDEx2",
                },
                self.comparative_dir / "10a_DEG_KO_ALDEx2" / "10_ALDEx2_KO_all_comparisons.tsv",
                "ALDEx2 KO comparative analysis",
            ),
            (
                "compare_ko_deseq2",
                self.comparative_dir / "10b_DEG_KO_DESeq2",
                {
                    "counts_file": ko_matrix,
                    "meta_file": meta_file,
                    "out_dir": self.comparative_dir / "10b_DEG_KO_DESeq2",
                    "ref_level": ref_level,
                },
                self.comparative_dir / "10b_DEG_KO_DESeq2" / "DESeq2_KO_all_comparisons.tsv",
                "DESeq2 KO comparative analysis",
            ),
            (
                "compare_pathway_aldex2",
                self.comparative_dir / "10a_DEG_KEGG_Pathway_ALDEx2",
                {
                    "counts_file": pathway_matrix,
                    "group_file": meta_file,
                    "out_dir": self.comparative_dir / "10a_DEG_KEGG_Pathway_ALDEx2",
                },
                self.comparative_dir / "10a_DEG_KEGG_Pathway_ALDEx2" / "ALDEx2_pairwise_all_comparisons.tsv",
                "ALDEx2 KEGG Pathway comparative analysis",
            ),
            (
                "compare_pathway_deseq2",
                self.comparative_dir / "10b_DEG_KEGG_Pathway_DESeq2",
                {
                    "counts_file": pathway_matrix,
                    "meta_file": meta_file,
                    "out_dir": self.comparative_dir / "10b_DEG_KEGG_Pathway_DESeq2",
                    "ref_level": ref_level,
                },
                self.comparative_dir / "10b_DEG_KEGG_Pathway_DESeq2" / "DESeq2_KEGG_Pathway_all_comparisons.tsv",
                "DESeq2 KEGG Pathway comparative analysis",
            ),
        ]

        for script_key, out_dir, assignments, expected, label in jobs:
            script = Path(SOURCE_SCRIPTS[script_key])
            if not script.exists():
                raise RuntimeError(f"Comparative source script missing: {script}")
            mkdir(out_dir)
            runnable = self.materialize_r_script(script, assignments, script_key)
            self.run_r_script(runnable, label)
            if not expected.exists():
                raise RuntimeError(f"{label} finished but expected output is missing: {expected}")

    def rnaseq_overview(self) -> None:
        script = Path(SOURCE_SCRIPTS["rnaseq_overview"])
        if not script.exists():
            raise RuntimeError(f"RNA-seq overview source script missing: {script}")
        configured_count_file = str(self.config.get("rnaseq_count_matrix") or "").strip()
        count_file = Path(configured_count_file or str(self.aggregate_counts_dir / "09b_COUNTS_KO_matrix.tsv")).expanduser()
        meta_file = self.outdir / "00_metadata.tsv"
        if not count_file.exists():
            raise RuntimeError(f"RNA-seq overview count matrix missing: {count_file}")
        if not meta_file.exists():
            raise RuntimeError(f"RNA-seq overview metadata missing: {meta_file}")
        header = read_table_header(count_file)
        if len(header) < 3:
            raise RuntimeError("RNA-seq overview count matrix must have a feature column and at least two sample columns")
        matrix_samples = [sample.strip() for sample in header[1:] if sample.strip()]
        duplicated_matrix_samples = duplicate_values(matrix_samples)
        if duplicated_matrix_samples:
            raise RuntimeError("RNA-seq overview count matrix has duplicated sample columns: " + ", ".join(duplicated_matrix_samples))
        metadata_samples = read_metadata_sample_names(meta_file)
        duplicated_metadata_samples = duplicate_values(metadata_samples)
        if duplicated_metadata_samples:
            raise RuntimeError("RNA-seq overview metadata has duplicated samples: " + ", ".join(duplicated_metadata_samples))
        missing_in_meta = sorted(set(matrix_samples) - set(metadata_samples))
        missing_in_counts = sorted(set(metadata_samples) - set(matrix_samples))
        if missing_in_meta:
            raise RuntimeError("Samples missing from RNA-seq overview metadata: " + ", ".join(missing_in_meta))
        if missing_in_counts:
            raise RuntimeError("Samples missing from RNA-seq overview count matrix: " + ", ".join(missing_in_counts))
        heatmap_sizes = parse_positive_int_list(self.config.get("rnaseq_heatmap_top", "50,100"), [50, 100])
        mkdir(self.rnaseq_dir)
        runnable_script = self.materialize_r_script(
            script,
            {
                "count_file": count_file,
                "meta_file": meta_file,
                "outdir": self.rnaseq_dir,
            },
            "13_RNAseq_one_group_analysis",
            heatmap_sizes=heatmap_sizes,
        )
        self.run_r_script(
            runnable_script,
            "RNA-seq one-group R overview",
        )
        expected_heatmaps = [self.rnaseq_dir / "figures" / f"Heatmap_vst_top{size}.pdf" for size in heatmap_sizes]
        expected = [
            self.rnaseq_dir / "tables" / "PCA_coordinates.tsv",
            self.rnaseq_dir / "tables" / "Gene_summary_mean_sd_cv.tsv",
            self.rnaseq_dir / "figures" / "PCA_vst.pdf",
        ] + expected_heatmaps
        missing_outputs = [str(path) for path in expected if not path.exists()]
        if missing_outputs:
            raise RuntimeError("RNA-seq overview finished but expected outputs are missing: " + ", ".join(missing_outputs))
        expected_heatmap_names = {path.name for path in expected_heatmaps}
        for heatmap in (self.rnaseq_dir / "figures").glob("Heatmap_vst_top*.pdf"):
            if heatmap.name not in expected_heatmap_names:
                heatmap.unlink()
        stale_readme = self.rnaseq_dir / "tables" / "README.txt"
        if stale_readme.exists() and "did not finish" in stale_readme.read_text(errors="replace"):
            stale_readme.unlink()

    def amr_ko(self) -> None:
        ko_matrix = self.aggregate_counts_dir / "09b_COUNTS_KO_matrix.tsv"
        if not ko_matrix.exists():
            raise RuntimeError(f"KO count matrix missing: {ko_matrix}")
        mkdir(self.amr_ko_dir / "tables")
        mkdir(self.amr_ko_dir / "figures")
        script = Path(SOURCE_SCRIPTS["amr_ko"])
        runnable_script = self.materialize_r_script(
            script,
            {
                "count_file": ko_matrix,
                "meta_file": self.outdir / "00_metadata.tsv",
                "outdir": self.amr_ko_dir,
            },
            "14_AMR_KO",
        )
        self.run_r_script(
            runnable_script,
            "KO AMR R analysis",
            allow_fail=True,
        )

        amr_kos = [
            "K18138",
            "K18139",
            "K18140",
            "K08204",
            "K11744",
            "K00817",
            "K01467",
            "K03883",
            "K18220",
            "K02542",
            "K07690",
            "K07787",
            "K03466",
        ]
        _first, samples, matrix = read_feature_matrix(ko_matrix)
        present = {ko: matrix[ko] for ko in amr_kos if ko in matrix}
        with (self.amr_ko_dir / "tables" / "AMR_KO_summary.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["KO", "sum_counts", "mean_counts", "samples_present"])
            for ko, values in sorted(present.items(), key=lambda item: sum(item[1]), reverse=True):
                writer.writerow([ko, int(sum(values)), f"{sum(values) / len(values):.3f}", sum(1 for value in values if value > 0)])
        activity = [0.0] * len(samples)
        for values in present.values():
            for idx, value in enumerate(values):
                activity[idx] += value
        with (self.amr_ko_dir / "tables" / "AMR_activity_index.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "AMR_activity"])
            for sample, value in zip(samples, activity):
                writer.writerow([sample, int(round(value))])
        if not present:
            (self.amr_ko_dir / "tables" / "NO_AMR_KO_DETECTED.txt").write_text("No configured AMR KO IDs were detected in the KO matrix.\n", encoding="utf-8")

    def card_rgi(self) -> None:
        tpm_file = self.salmon_matrix_dir / "salmon_TPM_all.tsv"
        if not tpm_file.exists():
            raise RuntimeError(f"TPM matrix missing: {tpm_file}")
        _first, samples, tpm_matrix = read_feature_matrix(tpm_file)
        threshold = float(self.config.get("card_tpm_threshold", 1.0))
        threads = int(self.config.get("threads", 16))
        failures = 0
        for sample_idx, sample in enumerate(samples):
            sample_out = self.card_dir / sample
            mkdir(sample_out)
            ids_p = sample_out / f"{sample}_TPM_gt{threshold}.p_ids.txt"
            contig_ids = sample_out / f"{sample}_TPM_gt{threshold}_contigs.txt"
            with ids_p.open("w") as handle:
                for gene, values in tpm_matrix.items():
                    if values[sample_idx] > threshold:
                        handle.write(gene + "\n")
            ids = ids_p.read_text().splitlines()
            if not ids:
                self.log(f"SKIP CARD {sample}: no transcripts above TPM>{threshold}")
                continue
            sample_prefix = f"{sample}__"
            with contig_ids.open("w") as handle:
                for contig in sorted({
                    re.sub(r"\.p[0-9]+$", "", gene[len(sample_prefix):] if gene.startswith(sample_prefix) else gene)
                    for gene in ids
                }):
                    handle.write(contig + "\n")
            asm_fasta = self.assembly_dir / sample / "final.contigs.fa"
            if not asm_fasta.exists():
                self.log(f"SKIP CARD {sample}: assembly FASTA missing: {asm_fasta}")
                continue
            expressed_fasta = sample_out / f"{sample}_TPM_gt{threshold}.fasta"
            script_extract = f"seqtk subseq {q(asm_fasta)} {q(contig_ids)} > {q(expressed_fasta)}"
            rc = self.run_bash("card", script_extract, f"CARD seqtk subseq {sample}", allow_fail=True)
            if rc != 0 or not expressed_fasta.exists() or expressed_fasta.stat().st_size == 0:
                failures += 1
                continue
            script_td = (
                f"TransDecoder.LongOrfs -t {q(expressed_fasta.name)}\n"
                f"TransDecoder.Predict -t {q(expressed_fasta.name)} --no_refine_starts"
            )
            rc = self.run_bash("meta", script_td, f"CARD TransDecoder {sample}", cwd=sample_out, allow_fail=True)
            pep = sample_out / f"{expressed_fasta.name}.transdecoder.pep"
            if rc != 0 or not pep.exists():
                failures += 1
                continue
            script_rgi = f"rgi -i {q(pep)} -o {q('card_amr_' + sample)} -t protein -n {threads} -e loose"
            rc = self.run_bash("card", script_rgi, f"RGI CARD {sample}", cwd=sample_out, allow_fail=True)
            if rc != 0 or not (sample_out / f"card_amr_{sample}.txt").exists():
                failures += 1
        if failures:
            raise RuntimeError(f"CARD/RGI failed for {failures} samples")

    def read_card_rows(self) -> list[dict[str, str]]:
        rows: list[dict[str, str]] = []
        for path in sorted(self.card_dir.glob("*/card_amr_*.txt")):
            sample = path.parent.name
            with path.open(newline="", errors="replace") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    cleaned = {
                        re.sub(r"[^a-z0-9]+", "_", key.strip().lower()).strip("_"): value
                        for key, value in row.items()
                        if key is not None
                    }
                    cleaned["sample"] = sample
                    cleaned["source_file"] = str(path)
                    rows.append(cleaned)
        return rows

    def card_summary(self) -> None:
        card_files = sorted(self.card_dir.glob("*/card_amr_*.txt"))
        if not card_files:
            raise RuntimeError(f"No CARD/RGI result files found under {self.card_dir}")
        rows = self.read_card_rows()
        cutoff_keep = {item.strip() for item in self.config.get("card_cutoffs", "Strict,Loose").split(",") if item.strip()}
        filtered = [row for row in rows if row.get("cut_off", "") in cutoff_keep]
        tables = self.card_summary_dir / "tables"
        mkdir(tables)
        if filtered:
            fieldnames = sorted({key for row in filtered for key in row.keys()})
        elif rows:
            fieldnames = sorted({key for row in rows for key in row.keys()})
        else:
            header = (card_files[0].read_text(errors="replace").splitlines() or [""])[0].split("\t")
            fieldnames = [re.sub(r"[^a-z0-9]+", "_", key.strip().lower()).strip("_") for key in header if key]
            fieldnames.extend(["sample", "source_file"])
        with (tables / "CARD_all_filtered.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(filtered)

        if not filtered:
            with (tables / "CARD_no_hits_or_no_hits_after_cutoff.tsv").open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["sample", "card_file", "raw_hits", "kept_hits", "note"])
                raw_counts: dict[str, int] = defaultdict(int)
                kept_counts: dict[str, int] = defaultdict(int)
                for row in rows:
                    raw_counts[row.get("sample", "")] += 1
                for row in filtered:
                    kept_counts[row.get("sample", "")] += 1
                for path in card_files:
                    sample = path.parent.name
                    note = "no CARD/RGI hits" if raw_counts[sample] == 0 else "hits filtered out by cutoff"
                    writer.writerow([sample, path, raw_counts[sample], kept_counts[sample], note])

        aro_counts: dict[tuple[str, str], int] = defaultdict(int)
        mech_counts: dict[tuple[str, str], int] = defaultdict(int)
        drug_counts: dict[tuple[str, str], int] = defaultdict(int)
        for row in filtered:
            sample = row.get("sample", "")
            aro = row.get("aro_name") or row.get("aro") or row.get("best_hit_aro")
            if aro:
                aro_counts[(sample, aro)] += 1
            category_value = row.get("best_hit_aro_category") or row.get("aro_category") or row.get("drug_class") or ""
            for raw_category in [part.strip() for part in category_value.split(";") if part.strip()]:
                lower = raw_category.lower()
                if "antibiotic inactivation" in lower:
                    mech = "inactivation"
                elif "efflux" in lower:
                    mech = "efflux"
                elif "target" in lower:
                    mech = "target"
                elif "reduced permeability" in lower:
                    mech = "permeability"
                else:
                    mech = "other"
                mech_counts[(sample, mech)] += 1
                if "resistance" in lower:
                    drug = re.sub(r"^determinant of\s+", "", raw_category, flags=re.I)
                    drug = re.sub(r"\s+resistance$", "", drug, flags=re.I).strip()
                    if drug:
                        drug_counts[(sample, drug)] += 1

        self.write_count_table(tables / "CARD_ARO_per_sample.tsv", ["sample", "aro_name", "n"], aro_counts)
        self.write_count_table(tables / "CARD_mechanisms_per_sample.tsv", ["sample", "mechanism", "n"], mech_counts)
        self.write_count_table(tables / "CARD_drug_classes_per_sample.tsv", ["sample", "drug_class", "n"], drug_counts)
        self.write_matrix_from_counts(tables / "CARD_mechanism_matrix.tsv", "mechanism", mech_counts)

    @staticmethod
    def write_count_table(path: Path, header: list[str], counts: dict[tuple[str, str], int]) -> None:
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(header)
            for (sample, feature), count in sorted(counts.items()):
                writer.writerow([sample, feature, count])

    @staticmethod
    def write_matrix_from_counts(path: Path, first_col: str, counts: dict[tuple[str, str], int]) -> None:
        samples = sorted({sample for sample, _feature in counts.keys()})
        features = sorted({feature for _sample, feature in counts.keys()})
        matrix = {feature: [float(counts.get((sample, feature), 0)) for sample in samples] for feature in features}
        write_feature_matrix(path, first_col, samples, matrix, integer=True)

    def card_ko_integration(self) -> None:
        card_table = self.card_summary_dir / "tables" / "CARD_all_filtered.tsv"
        ko_matrix = self.aggregate_counts_dir / "09b_COUNTS_KO_matrix.tsv"
        if not card_table.exists() or not ko_matrix.exists():
            raise RuntimeError("CARD summary table or KO matrix is missing")

        card_rows: list[dict[str, str]] = []
        with card_table.open(newline="", errors="replace") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            card_rows = list(reader)
        _first, ko_samples, ko_data = read_feature_matrix(ko_matrix)
        card_samples = sorted({row.get("sample", "") for row in card_rows if row.get("sample", "")})
        if not card_samples:
            card_samples = sorted({path.parent.name for path in self.card_dir.glob("*/card_amr_*.txt")})
        common = [sample for sample in ko_samples if sample in set(card_samples)]
        mkdir(self.integration_dir / "tables")
        if not common:
            (self.integration_dir / "tables" / "CARD_KO_integrated_metrics.tsv").write_text("sample\n", encoding="utf-8")
            self.log("No common samples between CARD and KO matrices")
            return

        aro_by_sample: dict[str, list[str]] = defaultdict(list)
        mech_by_sample: dict[str, list[str]] = defaultdict(list)
        drug_by_sample: dict[str, list[str]] = defaultdict(list)
        for row in card_rows:
            sample = row.get("sample", "")
            if sample not in common:
                continue
            aro = row.get("aro_name") or row.get("aro") or row.get("best_hit_aro")
            if aro:
                aro_by_sample[sample].append(aro)
            category_value = row.get("best_hit_aro_category") or row.get("aro_category") or row.get("drug_class") or ""
            for raw_category in [part.strip() for part in category_value.split(";") if part.strip()]:
                lower = raw_category.lower()
                if "antibiotic inactivation" in lower:
                    mech_by_sample[sample].append("inactivation")
                elif "efflux" in lower:
                    mech_by_sample[sample].append("efflux")
                elif "target" in lower:
                    mech_by_sample[sample].append("target")
                elif "reduced permeability" in lower:
                    mech_by_sample[sample].append("permeability")
                else:
                    mech_by_sample[sample].append("other")
                if "resistance" in lower:
                    drug = re.sub(r"^determinant of\s+", "", raw_category, flags=re.I)
                    drug = re.sub(r"\s+resistance$", "", drug, flags=re.I).strip()
                    if drug:
                        drug_by_sample[sample].append(drug)

        sample_to_idx = {sample: idx for idx, sample in enumerate(ko_samples)}
        metrics: dict[str, dict[str, float]] = {}
        for sample in common:
            idx = sample_to_idx[sample]
            metrics[sample] = {
                "amr_total_aro_hits": float(len(aro_by_sample[sample])),
                "amr_total_drug_hits": float(len(drug_by_sample[sample])),
                "amr_total_mech_hits": float(len(mech_by_sample[sample])),
                "amr_n_unique_aro": float(len(set(aro_by_sample[sample]))),
                "amr_n_unique_drugs": float(len(set(drug_by_sample[sample]))),
                "amr_n_unique_mech": float(len(set(mech_by_sample[sample]))),
                "ko_total_abund": float(sum(values[idx] for values in ko_data.values())),
                "ko_n_present": float(sum(1 for values in ko_data.values() if values[idx] > 0)),
            }

        metric_names = [
            "amr_total_aro_hits",
            "amr_total_drug_hits",
            "amr_total_mech_hits",
            "amr_n_unique_aro",
            "amr_n_unique_drugs",
            "amr_n_unique_mech",
            "ko_total_abund",
            "ko_n_present",
        ]
        with (self.integration_dir / "tables" / "CARD_KO_integrated_metrics.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample"] + metric_names)
            for sample in common:
                writer.writerow([sample] + [f"{metrics[sample][name]:.6g}" for name in metric_names])
        with (self.integration_dir / "tables" / "Spearman_correlations_AMR_KO_metrics.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["feature_1", "feature_2", "spearman_r"])
            vectors = {name: [metrics[sample][name] for sample in common] for name in metric_names}
            for name1 in metric_names:
                for name2 in metric_names:
                    writer.writerow([name1, name2, f"{spearman(vectors[name1], vectors[name2]):.6g}"])

    def final_multiqc(self) -> None:
        out = self.outdir / "99_multiqc_all"
        mkdir(out)
        inputs = [
            self.trim_dir / "reports",
            self.trim_dir / "multiqc",
            self.quast_dir,
            self.salmon_dir,
            self.card_dir,
        ]
        existing = [path for path in inputs if path.exists()]
        if not existing:
            self.log("SKIP MultiQC: no input directories")
            return
        script = "multiqc --force -o {out} -n multiqc_all.html {inputs}".format(
            out=q(out),
            inputs=" ".join(q(path) for path in existing),
        )
        self.run_bash("qc", script, "Final MultiQC", allow_fail=True)

    def report(self) -> None:
        rows = []
        suffixes = {".html", ".tsv", ".csv", ".txt", ".log", ".json", ".pdf", ".png"}
        for path in sorted(self.outdir.rglob("*")):
            if path.is_file() and path.suffix.lower() in suffixes:
                rel = path.relative_to(self.outdir)
                rows.append(
                    "<tr>"
                    f"<td>{html.escape(str(rel))}</td>"
                    f"<td>{path.stat().st_size}</td>"
                    "</tr>"
                )

        step_rows = []
        statuses = read_pipeline_step_statuses(self.outdir, list(self.steps))
        for status in statuses:
            step_rows.append(
                "<tr>"
                f"<td>{html.escape(status['step'])}</td>"
                f"<td>{html.escape(status['label'])}</td>"
                f"<td>{html.escape(status['status'])}</td>"
                f"<td>{html.escape(status.get('time', ''))}</td>"
                "</tr>"
            )

        source_rows = []
        for name, path in SOURCE_SCRIPTS.items():
            source_rows.append(f"<tr><td>{html.escape(name)}</td><td>{html.escape(path)}</td></tr>")

        cfg = html.escape(json.dumps(self.config, indent=2, ensure_ascii=True))
        report = self.outdir / "index.html"
        report.write_text(
            f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>RNA-seq AMR pipeline report</title>
<style>
body {{ font-family: system-ui, sans-serif; margin: 32px; color: #17202a; background: #f5f7fb; }}
h1 {{ margin-top: 0; }}
table {{ border-collapse: collapse; width: 100%; background: white; }}
td, th {{ border-bottom: 1px solid #d8dee9; padding: 7px 9px; text-align: left; vertical-align: top; }}
pre {{ background: #101820; color: #d6deeb; padding: 12px; border-radius: 8px; overflow: auto; }}
.grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 12px; margin-bottom: 20px; }}
.box {{ background: white; border: 1px solid #d8dee9; border-radius: 8px; padding: 12px; }}
</style>
</head>
<body>
<h1>RNA-seq AMR pipeline report</h1>
<div class="grid">
<div class="box"><b>Output</b><br>{html.escape(str(self.outdir))}</div>
<div class="box"><b>Samples</b><br>{len(self.samples)}</div>
<div class="box"><b>Finished</b><br>{html.escape(time.strftime("%Y-%m-%d %H:%M:%S"))}</div>
</div>
<h2>Steps</h2>
<table><thead><tr><th>Step</th><th>Label</th><th>Status</th><th>Time</th></tr></thead><tbody>{''.join(step_rows)}</tbody></table>
<h2>Configuration</h2>
<pre>{cfg}</pre>
<h2>Source scripts used as reference</h2>
<table><thead><tr><th>Name</th><th>Path</th></tr></thead><tbody>{''.join(source_rows)}</tbody></table>
<h2>Output files</h2>
<table><thead><tr><th>Path</th><th>Bytes</th></tr></thead><tbody>{''.join(rows)}</tbody></table>
</body>
</html>
""",
            encoding="utf-8",
        )
        self.log(f"Report written: {report}")


HTML_PAGE = """<!doctype html>
<html lang="pl">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>RNA-seq AMR Pipeline</title>
<style>
:root {
  color-scheme: light;
  --bg: #eef1f5;
  --panel: #ffffff;
  --panel-soft: #f8fafc;
  --line: #d5dde7;
  --text: #111827;
  --muted: #64748b;
  --blue: #1f5fbf;
  --blue-soft: #e8f1ff;
  --green: #157347;
  --green-soft: #e8f6ee;
  --red: #b42318;
  --red-soft: #fdebea;
  --amber: #996515;
  --amber-soft: #fff4d6;
  --shadow: 0 10px 28px rgba(15, 23, 42, .08);
}
* { box-sizing: border-box; }
body {
  margin: 0;
  min-width: 320px;
  background: var(--bg);
  color: var(--text);
  font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  font-size: 14px;
  letter-spacing: 0;
}
.topbar {
  min-height: 72px;
  padding: 14px 22px;
  background: var(--panel);
  border-bottom: 1px solid var(--line);
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 16px;
}
.brand h1 { margin: 0; font-size: 22px; line-height: 1.2; letter-spacing: 0; }
.brand-path { margin-top: 4px; color: var(--muted); font-size: 13px; }
.top-actions { display: flex; align-items: center; gap: 10px; flex-wrap: wrap; justify-content: flex-end; }
.badge {
  display: inline-flex;
  align-items: center;
  justify-content: center;
  min-height: 26px;
  padding: 3px 10px;
  border-radius: 999px;
  font-size: 12px;
  font-weight: 800;
  border: 1px solid transparent;
  background: #e5e7eb;
  color: #374151;
  white-space: nowrap;
}
.badge.running { background: var(--blue-soft); color: #174ea6; border-color: #bfd5ff; }
.badge.done { background: var(--green-soft); color: var(--green); border-color: #bfe6cc; }
.badge.failed { background: var(--red-soft); color: var(--red); border-color: #fac5bf; }
.layout {
  display: grid;
  grid-template-columns: minmax(460px, 600px) minmax(520px, 1fr);
  gap: 18px;
  padding: 18px;
  align-items: start;
}
.panel {
  background: var(--panel);
  border: 1px solid var(--line);
  border-radius: 8px;
  box-shadow: var(--shadow);
  min-width: 0;
}
.panel-head {
  min-height: 54px;
  padding: 14px 16px;
  border-bottom: 1px solid var(--line);
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 12px;
}
.panel-head h2 { margin: 0; font-size: 16px; line-height: 1.2; letter-spacing: 0; }
.panel-body { padding: 16px; }
.config-panel { max-height: calc(100vh - 108px); overflow: auto; }
.monitor-panel { min-height: calc(100vh - 108px); display: grid; grid-template-rows: auto auto auto 1fr; }
.form-section {
  padding: 0 0 16px;
  margin: 0 0 16px;
  border: 0;
  border-bottom: 1px solid var(--line);
}
.form-section:last-of-type { border-bottom: 0; margin-bottom: 0; padding-bottom: 0; }
legend {
  padding: 0;
  margin: 0 0 10px;
  font-size: 12px;
  text-transform: uppercase;
  font-weight: 850;
  color: #344054;
}
.field-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 12px; }
.field-grid.one { grid-template-columns: 1fr; }
label.field { display: grid; gap: 5px; margin: 0; font-weight: 700; color: #263244; font-size: 12px; }
input[type=text], input[type=number], select {
  width: 100%;
  min-height: 38px;
  padding: 8px 10px;
  border: 1px solid #bac5d3;
  border-radius: 6px;
  background: #fff;
  color: var(--text);
  font: inherit;
}
input:focus, select:focus { outline: 2px solid rgba(31, 95, 191, .18); border-color: var(--blue); }
.segmented { display: grid; grid-template-columns: 1fr 1fr; gap: 8px; }
.segmented label {
  display: flex;
  align-items: center;
  gap: 8px;
  min-height: 40px;
  padding: 9px 10px;
  border: 1px solid var(--line);
  border-radius: 6px;
  background: var(--panel-soft);
  font-weight: 750;
}
.test-panel {
  margin-top: 12px;
  padding: 12px;
  border: 1px solid #bfd5ff;
  border-left: 4px solid var(--blue);
  border-radius: 6px;
  background: #f7fbff;
}
.checks {
  display: grid;
  grid-template-columns: 1fr 1fr;
  gap: 7px 10px;
}
.checks label {
  display: flex;
  align-items: center;
  gap: 8px;
  min-height: 32px;
  padding: 6px 8px;
  border: 1px solid var(--line);
  border-radius: 6px;
  background: var(--panel-soft);
  font-size: 12px;
  font-weight: 650;
}
input[type=checkbox], input[type=radio] { accent-color: var(--blue); }
.recursive-row {
  margin-top: 12px;
  display: flex;
  align-items: center;
  gap: 8px;
  font-weight: 750;
  color: #263244;
}
.actions {
  position: sticky;
  bottom: 0;
  display: flex;
  gap: 10px;
  padding: 12px 0 0;
  margin-top: 16px;
  background: linear-gradient(to bottom, rgba(255,255,255,0), var(--panel) 24%);
}
button {
  min-height: 40px;
  padding: 9px 14px;
  border: 1px solid transparent;
  border-radius: 6px;
  background: var(--blue);
  color: white;
  font-weight: 850;
  cursor: pointer;
}
button.secondary { background: #f1f5f9; color: #162033; border-color: #cbd5e1; }
button:hover { filter: brightness(.98); }
.status-grid {
  display: grid;
  grid-template-columns: repeat(3, minmax(0, 1fr));
  gap: 10px;
  padding: 16px;
  border-bottom: 1px solid var(--line);
}
.metric {
  min-height: 70px;
  padding: 10px 12px;
  border: 1px solid var(--line);
  border-radius: 8px;
  background: var(--panel-soft);
}
.metric-label { color: var(--muted); font-size: 11px; text-transform: uppercase; font-weight: 850; }
.metric-value { margin-top: 6px; font-size: 14px; font-weight: 800; overflow-wrap: anywhere; }
.preview-box {
  margin: 0 16px 14px;
  min-height: 38px;
  padding: 10px 12px;
  border: 1px solid var(--line);
  border-radius: 8px;
  background: #fbfdff;
  color: var(--muted);
  overflow-wrap: anywhere;
}
.stages {
  padding: 0 16px 14px;
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(235px, 1fr));
  gap: 8px;
}
.stage {
  display: grid;
  grid-template-columns: 24px 1fr;
  gap: 8px;
  align-items: start;
  min-height: 56px;
  padding: 8px 9px;
  border: 1px solid var(--line);
  border-radius: 8px;
  background: var(--panel-soft);
}
.stage .icon {
  width: 22px;
  height: 22px;
  border-radius: 50%;
  display: inline-grid;
  place-items: center;
  font-size: 12px;
  font-weight: 900;
  background: #e2e8f0;
  color: #475569;
}
.stage.done .icon { background: var(--green-soft); color: var(--green); }
.stage.running .icon { background: var(--blue-soft); color: var(--blue); }
.stage.failed .icon { background: var(--red-soft); color: var(--red); }
.stage.pending .icon { background: var(--amber-soft); color: var(--amber); }
.stage.inactive { opacity: .42; }
.stage-title { font-size: 12px; font-weight: 850; line-height: 1.25; }
.stage-meta { margin-top: 3px; font-size: 11px; color: var(--muted); overflow-wrap: anywhere; }
.log-wrap { padding: 0 16px 16px; min-height: 0; }
pre {
  margin: 0;
  min-height: 420px;
  max-height: calc(100vh - 430px);
  overflow: auto;
  white-space: pre-wrap;
  padding: 12px;
  border-radius: 8px;
  background: #0f172a;
  color: #dbeafe;
  font: 12px/1.45 ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
}
.small { color: var(--muted); font-size: 12px; }
@media (max-width: 1180px) {
  .layout { grid-template-columns: 1fr; }
  .config-panel, .monitor-panel { max-height: none; min-height: 0; }
  pre { max-height: 520px; }
}
@media (max-width: 720px) {
  .topbar { align-items: flex-start; flex-direction: column; }
  .layout { padding: 10px; gap: 10px; }
  .field-grid, .checks, .segmented, .status-grid { grid-template-columns: 1fr; }
  .panel-body, .panel-head, .status-grid, .log-wrap { padding-left: 12px; padding-right: 12px; }
  .preview-box { margin-left: 12px; margin-right: 12px; }
  .stages { padding-left: 12px; padding-right: 12px; grid-template-columns: 1fr; }
}
</style>
</head>
<body>
<header class="topbar">
  <div class="brand">
    <h1>RNA-seq AMR Pipeline</h1>
    <div class="brand-path">FASTQ PE -> QC -> assembly -> quantification -> annotation -> AMR</div>
  </div>
  <div class="top-actions">
    <span id="statusBadge" class="badge">idle</span>
    <span class="badge">local</span>
  </div>
</header>
<main class="layout">
  <section class="panel config-panel">
    <div class="panel-head">
      <h2>Konfiguracja runu</h2>
      <span id="selectedCount" class="badge">0 steps</span>
    </div>
    <form id="runForm" class="panel-body">
      <fieldset class="form-section">
        <legend>Sciezki</legend>
        <div class="field-grid one">
          <label class="field">Input FASTQ directory
            <input type="text" name="input_dir" value="./input_fastq">
          </label>
          <label class="field">Output directory
            <input type="text" name="output_dir" value="./runs/run1">
          </label>
          <label class="field">Metadata TSV/CSV
            <input type="text" name="metadata_file" value="">
          </label>
        </div>
      </fieldset>

      <fieldset class="form-section">
        <legend>RNA i porownania</legend>
        <div class="field-grid">
          <label class="field">RNA-seq count matrix
            <input type="text" name="rnaseq_count_matrix" value="">
          </label>
          <label class="field">Heatmap top genes
            <input type="text" name="rnaseq_heatmap_top" value="50,100">
          </label>
          <label class="field">Grouping column
            <input type="text" name="comparative_group_col" value="group_all">
          </label>
          <label class="field">Reference group
            <input type="text" name="comparative_ref_level" value="">
          </label>
          <label class="field">Rscript binary
            <input type="text" name="rscript_bin" value="/usr/bin/Rscript">
          </label>
        </div>
      </fieldset>

      <fieldset class="form-section">
        <legend>Run profile</legend>
        <div class="segmented">
          <label><input type="radio" name="run_profile" value="full" checked> Full run</label>
          <label><input type="radio" name="run_profile" value="test"> Test run</label>
        </div>
        <div id="testPanel" class="test-panel" hidden>
          <div class="field-grid">
            <label class="field">Test samples
              <input type="number" name="test_sample_count" value="2" min="1">
            </label>
            <label class="field">Reads per sample
              <input type="number" name="test_read_count" value="10000" min="1">
            </label>
            <label class="field">Random seed
              <input type="number" name="test_seed" value="7" min="1">
            </label>
          </div>
        </div>
      </fieldset>

      <fieldset class="form-section">
        <legend>Resources and databases</legend>
        <div class="field-grid">
          <label class="field">Threads
            <input type="number" name="threads" value="16" min="1">
          </label>
          <label class="field">Salmon threads
            <input type="number" name="salmon_threads" value="8" min="1">
          </label>
          <label class="field">Bowtie2 rRNA index
            <input type="text" name="rrna_index" value="">
          </label>
          <label class="field">eggNOG data directory
            <input type="text" name="eggnog_data_dir" value="">
          </label>
        </div>
      </fieldset>

      <fieldset class="form-section">
        <legend>Parametry</legend>
        <div class="field-grid">
          <label class="field">fastp strategy
            <select name="fastp_strategy">
              <option value="full">adapters + ends + quality</option>
              <option value="adapters">adapters only</option>
              <option value="ends">adapters + cut_front/cut_tail</option>
              <option value="quality">adapters + Phred/min length</option>
              <option value="custom">custom fastp options</option>
            </select>
          </label>
          <label class="field">QUAST mode
            <select name="quast_mode">
              <option value="quast">QUAST de novo report</option>
              <option value="metaquast">MetaQUAST no reference search</option>
            </select>
          </label>
          <label class="field">fastp Phred
            <input type="number" name="fastp_phred" value="20" min="1">
          </label>
          <label class="field">fastp min length
            <input type="number" name="fastp_min_len" value="50" min="1">
          </label>
          <label class="field">Custom fastp options
            <input type="text" name="fastp_extra" value="--detect_adapter_for_pe">
          </label>
          <label class="field">MEGAHIT k-list
            <input type="text" name="megahit_k_list" value="21,41,61,81,101,121">
          </label>
          <label class="field">MEGAHIT min contig length
            <input type="number" name="megahit_min_contig_len" value="300" min="100">
          </label>
          <label class="field">TPM threshold
            <input type="number" name="tpm_threshold" value="1" min="0" step="0.1">
          </label>
          <label class="field">Min samples for TPM filter
            <input type="number" name="tpm_min_samples" value="1" min="1">
          </label>
          <label class="field">CARD TPM threshold
            <input type="number" name="card_tpm_threshold" value="1" min="0" step="0.1">
          </label>
          <label class="field">CARD cutoffs
            <input type="text" name="card_cutoffs" value="Strict,Loose">
          </label>
        </div>
      </fieldset>

      <fieldset class="form-section">
        <legend>Execution scope</legend>
        <div class="field-grid">
          <label class="field">Start from step
            <select name="start_from">
              <option value="beginning">beginning</option>
              <option value="estimate_rrna">estimate rRNA</option>
              <option value="trim_qc">fastp</option>
              <option value="remove_rrna">remove rRNA</option>
              <option value="megahit">MEGAHIT</option>
              <option value="transdecoder">TransDecoder</option>
              <option value="salmon">Salmon</option>
              <option value="eggnog">eggNOG</option>
              <option value="aggregate">aggregate</option>
              <option value="kegg_pathways">KEGG pathways</option>
              <option value="comparative">comparative analysis</option>
              <option value="rnaseq_overview">RNA-seq overview</option>
              <option value="amr_ko">AMR KO</option>
              <option value="card_rgi">CARD/RGI</option>
            </select>
          </label>
          <label class="field">Run mode
            <select name="run_mode">
              <option value="resume">resume completed steps</option>
              <option value="force">force rerun selected steps</option>
            </select>
          </label>
        </div>
      </fieldset>

      <fieldset class="form-section">
        <legend>Steps</legend>
        <div class="checks">
          <label><input type="checkbox" name="steps" value="estimate_rrna"> estimate rRNA content</label>
          <label><input type="checkbox" name="steps" value="trim_qc" checked> fastp trim + MultiQC</label>
          <label><input type="checkbox" name="steps" value="remove_rrna" checked> remove rRNA Bowtie2</label>
          <label><input type="checkbox" name="steps" value="rrna_stats" checked> rRNA removal stats</label>
          <label><input type="checkbox" name="steps" value="megahit" checked> MEGAHIT assembly</label>
          <label><input type="checkbox" name="steps" value="collect_contigs" checked> collect contigs</label>
          <label><input type="checkbox" name="steps" value="quast" checked> QUAST / MetaQUAST report</label>
          <label><input type="checkbox" name="steps" value="transdecoder" checked> TransDecoder CDS/PEP</label>
          <label><input type="checkbox" name="steps" value="collect_cds" checked> collect CDS/proteins</label>
          <label><input type="checkbox" name="steps" value="cds_stats" checked> CDS stats</label>
          <label><input type="checkbox" name="steps" value="salmon" checked> Salmon quant</label>
          <label><input type="checkbox" name="steps" value="salmon_matrix" checked> Salmon matrices</label>
          <label><input type="checkbox" name="steps" value="eggnog" checked> eggNOG mapper</label>
          <label><input type="checkbox" name="steps" value="aggregate" checked> KO/EC/PFAM aggregation</label>
          <label><input type="checkbox" name="steps" value="kegg_pathways" checked> KO to KEGG pathways</label>
          <label><input type="checkbox" name="steps" value="comparative" checked> comparative KO/Pathway</label>
          <label><input type="checkbox" name="steps" value="rnaseq_overview" checked> RNA-seq overview R plots</label>
          <label><input type="checkbox" name="steps" value="amr_ko" checked> AMR from KO matrix</label>
          <label><input type="checkbox" name="steps" value="card_rgi" checked> CARD/RGI</label>
          <label><input type="checkbox" name="steps" value="card_summary" checked> CARD summary</label>
          <label><input type="checkbox" name="steps" value="card_ko_integration" checked> CARD x KO integration</label>
          <label><input type="checkbox" name="steps" value="multiqc" checked> final MultiQC</label>
          <label><input type="checkbox" name="steps" value="summary" checked> final HTML index</label>
        </div>
        <label class="recursive-row"><input type="checkbox" name="recursive" checked> recursive file search</label>
      </fieldset>

      <div class="actions">
        <button type="button" id="previewBtn" class="secondary">Preview samples</button>
        <button type="submit">Start pipeline</button>
      </div>
    </form>
  </section>

  <section class="panel monitor-panel">
    <div class="panel-head">
      <h2>Monitoring</h2>
      <span class="small" id="jobLabel">Idle</span>
    </div>
    <div class="status-grid">
      <div class="metric">
        <div class="metric-label">Status</div>
        <div class="metric-value" id="status">Idle</div>
      </div>
      <div class="metric">
        <div class="metric-label">Output</div>
        <div class="metric-value" id="outputMetric">-</div>
      </div>
      <div class="metric">
        <div class="metric-label">Config</div>
        <div class="metric-value" id="configMetric">-</div>
      </div>
    </div>
    <div id="preview" class="preview-box">No sample preview yet.</div>
    <div id="stages" class="stages"></div>
    <div class="log-wrap"><pre id="log"></pre></div>
  </section>
</main>
<script>
const form = document.getElementById('runForm');
const logEl = document.getElementById('log');
const statusEl = document.getElementById('status');
const statusBadge = document.getElementById('statusBadge');
const previewEl = document.getElementById('preview');
const previewBtn = document.getElementById('previewBtn');
const stagesEl = document.getElementById('stages');
const testPanel = document.getElementById('testPanel');
const selectedCount = document.getElementById('selectedCount');
const jobLabel = document.getElementById('jobLabel');
const outputMetric = document.getElementById('outputMetric');
const configMetric = document.getElementById('configMetric');
let jobId = null;
const icons = { done: '*', running: '>', failed: '!', pending: '.', inactive: '-' };
function escapeHtml(value) {
  return String(value || '').replace(/[&<>"']/g, ch => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[ch]));
}
function updateMode() {
  const profile = new FormData(form).get('run_profile');
  testPanel.hidden = profile !== 'test';
}
function updateSelectedCount() {
  const n = form.querySelectorAll('input[name=steps]:checked').length;
  selectedCount.textContent = n + ' steps';
}
form.querySelectorAll('input[name=run_profile]').forEach(el => el.addEventListener('change', updateMode));
form.querySelectorAll('input[name=steps]').forEach(el => el.addEventListener('change', updateSelectedCount));
updateMode();
updateSelectedCount();
function renderStages(steps) {
  if (!steps || !steps.length) {
    stagesEl.innerHTML = '';
    return;
  }
  stagesEl.innerHTML = steps.map(step => {
    const status = step.status || 'pending';
    const meta = step.error ? step.error : (step.time || status);
    return `<div class="stage ${escapeHtml(status)}">
      <span class="icon">${icons[status] || '.'}</span>
      <div><div class="stage-title">${escapeHtml(step.label)}</div><div class="stage-meta">${escapeHtml(meta)}</div></div>
    </div>`;
  }).join('');
}
function setBadge(status) {
  statusBadge.className = 'badge ' + (status || '');
  statusBadge.textContent = status || 'idle';
}
previewBtn.addEventListener('click', async () => {
  const data = new FormData(form);
  const params = new URLSearchParams();
  params.set('input_dir', data.get('input_dir'));
  params.set('recursive', data.get('recursive') ? 'true' : 'false');
  previewEl.textContent = 'Scanning input directory...';
  const res = await fetch('/preview?' + params.toString());
  const js = await res.json();
  if (!js.ok) {
    previewEl.textContent = 'Preview error: ' + js.error;
    return;
  }
  const names = js.samples.map(s => s.sample).join(', ');
  previewEl.textContent = 'Detected paired samples: ' + js.count + (names ? ' | ' + names : '');
});
form.addEventListener('submit', async (event) => {
  event.preventDefault();
  const data = new FormData(form);
  const res = await fetch('/run', { method: 'POST', body: new URLSearchParams(data) });
  const js = await res.json();
  if (!res.ok || js.error) {
    statusEl.textContent = 'Start error';
    outputMetric.textContent = js.error || res.statusText;
    configMetric.textContent = '-';
    setBadge('failed');
    return;
  }
  jobId = js.job_id;
  jobLabel.textContent = 'Job ' + jobId;
  statusEl.textContent = 'running';
  outputMetric.textContent = js.output_dir || '-';
  configMetric.textContent = js.config || '-';
  setBadge('running');
  poll();
});
async function poll() {
  if (!jobId) return;
  const res = await fetch('/status?job_id=' + encodeURIComponent(jobId));
  const js = await res.json();
  statusEl.textContent = js.status || 'unknown';
  outputMetric.textContent = js.output_dir || '-';
  configMetric.textContent = js.config || '-';
  setBadge(js.status);
  renderStages(js.steps || []);
  logEl.textContent = js.log || '';
  setTimeout(poll, js.done ? 5000 : 2000);
}
</script>
</body>
</html>
"""


JOBS: dict[str, dict[str, Any]] = {}


def read_pipeline_step_statuses(outdir: Path, selected_steps: list[str]) -> list[dict[str, str]]:
    selected = set(selected_steps)
    state_dir = outdir / ".pipeline_state"
    statuses: list[dict[str, str]] = []
    for step, label in PIPELINE_STEPS:
        status = "inactive" if step not in selected else "pending"
        timestamp = ""
        error = ""
        state_json = state_dir / f"{step}.json"
        if state_json.exists():
            try:
                payload = json.loads(state_json.read_text(errors="replace"))
                timestamp = str(payload.get("time", ""))
                error = str(payload.get("error", ""))
            except json.JSONDecodeError:
                pass
        if (state_dir / f"{step}.failed").exists():
            status = "failed"
        elif (state_dir / f"{step}.running").exists():
            status = "running"
        elif (state_dir / f"{step}.done").exists():
            status = "done"
        statuses.append({"step": step, "label": label, "status": status, "time": timestamp, "error": error})
    return statuses


class Handler(BaseHTTPRequestHandler):
    def _send(self, status: int, body: bytes, content_type: str) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        if parsed.path == "/":
            self._send(200, HTML_PAGE.encode("utf-8"), "text/html; charset=utf-8")
            return
        if parsed.path == "/preview":
            qs = parse_qs(parsed.query)
            input_dir = Path(qs.get("input_dir", [""])[0]).expanduser()
            recursive = qs.get("recursive", ["true"])[0] == "true"
            try:
                samples = discover_samples(input_dir, recursive=recursive)
                payload = {
                    "ok": True,
                    "count": len(samples),
                    "samples": [
                        {"sample": s.sample, "r1": str(s.r1 or ""), "r2": str(s.r2 or "")}
                        for s in samples[:100]
                    ],
                }
            except Exception as exc:
                payload = {"ok": False, "error": str(exc), "count": 0, "samples": []}
            self._send(200, json.dumps(payload).encode("utf-8"), "application/json")
            return
        if parsed.path == "/status":
            qs = parse_qs(parsed.query)
            job_id = qs.get("job_id", [""])[0]
            job = JOBS.get(job_id)
            if not job:
                self._send(404, b'{"error":"unknown job"}', "application/json")
                return
            log_path = Path(job["log"])
            log = ""
            if log_path.exists():
                log = log_path.read_text(errors="replace")[-50000:]
            payload = {
                "job_id": job_id,
                "status": job["status"],
                "done": job["status"] in {"done", "failed"},
                "log": log,
                "config": job["config"],
                "output_dir": job["output_dir"],
                "steps": read_pipeline_step_statuses(Path(job["output_dir"]), job.get("steps", [])),
            }
            self._send(200, json.dumps(payload).encode("utf-8"), "application/json")
            return
        self._send(404, b"not found", "text/plain")

    def do_POST(self) -> None:
        parsed = urlparse(self.path)
        if parsed.path != "/run":
            self._send(404, b"not found", "text/plain")
            return
        length = int(self.headers.get("Content-Length", "0"))
        body = self.rfile.read(length).decode("utf-8")
        form = parse_qs(body)
        try:
            config = form_to_config(form)
        except Exception as exc:
            self._send(400, json.dumps({"error": str(exc)}).encode("utf-8"), "application/json")
            return
        outdir = Path(config["output_dir"]).expanduser().resolve()
        mkdir(outdir)
        config_path = outdir / "pipeline_config.json"
        config_path.write_text(json.dumps(config, indent=2), encoding="utf-8")
        job_id = str(int(time.time() * 1000))
        log_path = outdir / "pipeline.log"
        JOBS[job_id] = {
            "status": "running",
            "log": str(log_path),
            "config": str(config_path),
            "output_dir": str(outdir),
            "steps": config.get("steps", []),
        }

        def worker() -> None:
            cmd = [sys.executable, str(Path(__file__).resolve()), "run", "--config", str(config_path)]
            proc = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, text=True)
            JOBS[job_id]["status"] = "done" if proc.returncode == 0 else "failed"

        threading.Thread(target=worker, daemon=True).start()
        payload = {"job_id": job_id, "output_dir": str(outdir), "config": str(config_path)}
        self._send(200, json.dumps(payload).encode("utf-8"), "application/json")


def form_to_config(form: dict[str, list[str]]) -> dict[str, Any]:
    def one(key: str, default: str = "") -> str:
        return form.get(key, [default])[0]

    input_dir = one("input_dir").strip()
    output_dir = one("output_dir").strip()
    if not input_dir:
        raise ValueError("Input directory is empty")
    if not output_dir:
        raise ValueError("Output directory is empty")

    return {
        "input_dir": input_dir,
        "output_dir": output_dir,
        "metadata_file": one("metadata_file", "").strip(),
        "rnaseq_count_matrix": one("rnaseq_count_matrix", "").strip(),
        "rnaseq_heatmap_top": one("rnaseq_heatmap_top", "50,100").strip(),
        "comparative_group_col": one("comparative_group_col", "group_all").strip(),
        "comparative_ref_level": one("comparative_ref_level", "").strip(),
        "rscript_bin": one("rscript_bin", "/usr/bin/Rscript").strip(),
        "run_profile": one("run_profile", "full"),
        "test_sample_count": int(one("test_sample_count", "2") or "2"),
        "test_read_count": int(one("test_read_count", "10000") or "10000"),
        "test_seed": int(one("test_seed", "7") or "7"),
        "threads": int(one("threads", "16") or "16"),
        "salmon_threads": int(one("salmon_threads", "8") or "8"),
        "rrna_index": one("rrna_index", DEFAULT_DBS["rrna_bowtie2"]),
        "eggnog_data_dir": one("eggnog_data_dir", DEFAULT_DBS["eggnog_data"]),
        "fastp_strategy": one("fastp_strategy", "full"),
        "fastp_phred": int(one("fastp_phred", "20") or "20"),
        "fastp_min_len": int(one("fastp_min_len", "50") or "50"),
        "fastp_extra": one("fastp_extra", "--detect_adapter_for_pe"),
        "megahit_k_list": one("megahit_k_list", "21,41,61,81,101,121"),
        "megahit_min_contig_len": int(one("megahit_min_contig_len", "300") or "300"),
        "quast_mode": one("quast_mode", "quast"),
        "tpm_threshold": float(one("tpm_threshold", "1") or "1"),
        "tpm_min_samples": int(one("tpm_min_samples", "1") or "1"),
        "card_tpm_threshold": float(one("card_tpm_threshold", "1") or "1"),
        "card_cutoffs": one("card_cutoffs", "Strict,Loose"),
        "resume": one("run_mode", "resume") == "resume",
        "force": one("run_mode", "resume") == "force",
        "start_from": one("start_from", "beginning"),
        "recursive": "recursive" in form,
        "steps": form.get("steps", []),
    }


def serve(args: argparse.Namespace) -> None:
    server = ThreadingHTTPServer((args.host, args.port), Handler)
    print(f"Serving RNA-seq AMR pipeline app at http://{args.host}:{args.port}", flush=True)
    server.serve_forever()


def run_config(args: argparse.Namespace) -> None:
    config = json.loads(Path(args.config).read_text())
    try:
        Pipeline(config).run()
    except Exception:
        outdir = Path(config.get("output_dir", ".")).expanduser().resolve()
        mkdir(outdir)
        with (outdir / "pipeline.log").open("a") as handle:
            handle.write(traceback.format_exc())
        raise


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="RNA-seq AMR pipeline app")
    sub = parser.add_subparsers(dest="cmd", required=True)
    serve_p = sub.add_parser("serve", help="start local HTML app")
    serve_p.add_argument("--host", default="127.0.0.1")
    serve_p.add_argument("--port", type=int, default=8791)
    serve_p.set_defaults(func=serve)
    run_p = sub.add_parser("run", help="run pipeline from JSON config")
    run_p.add_argument("--config", required=True)
    run_p.set_defaults(func=run_config)
    args = parser.parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
