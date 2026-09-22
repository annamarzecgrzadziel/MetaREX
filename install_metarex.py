#!/usr/bin/env python3
"""Install and validate the local MetaREX runtime.

The installer deliberately uses the repository's explicit Conda lockfiles.
Large reference databases are never downloaded implicitly: pass their paths
with --rrna-index/--eggnog-data/--card-data, or provide explicit archive URLs
with the corresponding --*-url options.

Examples
--------
  python3 install_metarex.py install --prefix ~/metarex
  python3 install_metarex.py install --prefix ~/metarex \
      --rrna-index /path/to/SILVA_138_2_rRNA --eggnog-data /path/to/eggnog
  python3 install_metarex.py check --config ~/metarex/metarex_config.json
  python3 install_metarex.py versions --config ~/metarex/metarex_config.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import urllib.request
from urllib.parse import urlparse
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parent

LOCKS = {
    "qc": "qc-linux-64.txt",
    "meta": "meta-linux-64.txt",
    "quast": "quast-linux-64.txt",
    "bioinfo": "quast-linux-64.txt",
    "eggnog": "eggnog-linux-64.txt",
    "card": "card-linux-64.txt",
    "r": "meta-linux-64.txt",
}

TOOLS = {
    "qc": ("fastp", "multiqc", "bowtie2", "megahit"),
    "meta": ("TransDecoder.LongOrfs", "TransDecoder.Predict", "salmon"),
    "quast": ("quast.py", "metaquast.py"),
    "bioinfo": ("multiqc",),
    "eggnog": ("emapper.py",),
    "card": ("rgi", "seqtk"),
    "r": ("Rscript",),
}

# R packages used by the scripts in translate/.  The lockfiles contain the
# base R runtime, but an explicit lockfile does not guarantee that these
# script-level packages are present.  Install them after the environments are
# created and validate them with Rscript.
R_PACKAGES = {
    "DESeq2": "bioconductor-deseq2",
    "ALDEx2": "bioconductor-aldex2",
    "apeglm": "bioconductor-apeglm",
    "ggplot2": "r-ggplot2",
    "matrixStats": "r-matrixstats",
    "pheatmap": "r-pheatmap",
    "tidyverse": "r-tidyverse",
    "janitor": "r-janitor",
    "FactoMineR": "r-factominer",
    "factoextra": "r-factoextra",
    "KEGGREST": "bioconductor-keggrest",
}

PYTHON_PACKAGES = {
    "pandas": "pandas",
    "Bio.SeqIO": "biopython",
}

def log(message: str) -> None:
    print(f"[MetaREX installer] {message}", flush=True)


def run(command: list[str], *, check: bool = True, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    log("$ " + " ".join(str(x) for x in command))
    return subprocess.run(command, text=True, check=check, env=env)


def find_lock_dir() -> Path:
    """Find lockfiles without depending on a repository-specific directory name."""
    candidates = [ROOT / "environment_lockfiles"]
    candidates.extend(path for path in ROOT.rglob("environment_lockfiles") if path.is_dir())
    required = set(LOCKS.values())
    for candidate in candidates:
        if required.issubset({path.name for path in candidate.glob("*.txt")}):
            return candidate
    raise RuntimeError("Could not locate the five Linux environment lockfiles")


def find_engine(requested: str) -> str:
    if requested != "auto":
        path = shutil.which(requested)
        if not path:
            raise RuntimeError(f"Requested package manager not found: {requested}")
        return path
    for candidate in ("mamba", "conda", "micromamba"):
        path = shutil.which(candidate)
        if path:
            return path
    raise RuntimeError("No conda, mamba or micromamba executable was found in PATH")


def engine_create(engine: str, prefix: Path, lockfile: Path, package_cache: Path) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    if prefix.exists() and (prefix / "conda-meta").exists():
        log(f"Environment already exists: {prefix}")
        return
    prefix.parent.mkdir(parents=True, exist_ok=True)
    package_cache.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["CONDA_PKGS_DIRS"] = str(package_cache)
    env["MAMBA_ROOT_PREFIX"] = str(package_cache.parent / "mamba_root")
    run([engine, "create", "-y", "-p", str(prefix), "--file", str(lockfile)], env=env)


def install_r_packages(engine: str, env: Path, package_cache: Path) -> None:
    """Install the R packages required by MetaREX analysis scripts."""
    package_cache.mkdir(parents=True, exist_ok=True)
    env_vars = os.environ.copy()
    env_vars["CONDA_PKGS_DIRS"] = str(package_cache)
    env_vars["MAMBA_ROOT_PREFIX"] = str(package_cache.parent / "mamba_root")
    packages = list(R_PACKAGES.values()) + list(PYTHON_PACKAGES.values())
    log("Installing analysis packages: " + ", ".join(packages))
    run(
        [
            engine,
            "install",
            "-y",
            "-p",
            str(env),
            "-c",
            "conda-forge",
            "-c",
            "bioconda",
            *packages,
        ],
        env=env_vars,
    )


def run_in_env(engine: str, env: Path, command: list[str], check: bool = False) -> subprocess.CompletedProcess[str]:
    return subprocess.run([engine, "run", "-p", str(env), *command], text=True, capture_output=True, check=check)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download_and_extract(url: str, destination: Path) -> Path:
    destination.mkdir(parents=True, exist_ok=True)
    archive = destination / Path(urlparse(url).path).name
    if not archive.name or archive.name == "/":
        archive = destination / "metarex_database_download"
    log(f"Downloading {url} -> {archive}")
    urllib.request.urlretrieve(url, archive)
    extracted = destination / "extracted"
    extracted.mkdir(exist_ok=True)
    if zipfile.is_zipfile(archive):
        with zipfile.ZipFile(archive) as handle:
            handle.extractall(extracted)
    elif tarfile.is_tarfile(archive):
        with tarfile.open(archive) as handle:
            root = extracted.resolve()
            for member in handle.getmembers():
                target = (extracted / member.name).resolve()
                if target != root and root not in target.parents:
                    raise RuntimeError(f"Unsafe archive member: {member.name}")
            handle.extractall(extracted)
    else:
        log(f"Downloaded file is not a recognised archive; keeping it at {archive}")
        return archive
    return extracted


def copy_runtime(prefix: Path) -> Path:
    """Copy the application runtime, without analysis source scripts."""
    runtime = prefix / "share" / "metarex"
    runtime.mkdir(parents=True, exist_ok=True)
    for filename in ("rnaseq_amr_pipeline_app.py", "example_pipeline_config.json", "README.md"):
        source = ROOT / filename
        if source.exists():
            shutil.copy2(source, runtime / filename)
    launcher = prefix / "bin" / "metarex"
    launcher.parent.mkdir(parents=True, exist_ok=True)
    launcher.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        f"exec python3 {runtime / 'rnaseq_amr_pipeline_app.py'} \"$@\"\n",
        encoding="utf-8",
    )
    launcher.chmod(0o755)
    return runtime


def check_path(label: str, path: Path, errors: list[str], *, directory: bool = True) -> None:
    exists = path.is_dir() if directory else path.exists()
    if exists:
        log(f"OK {label}: {path}")
    else:
        errors.append(f"Missing {label}: {path}")


def check_rrna_index(path: Path, errors: list[str]) -> None:
    """Bowtie2 accepts an index prefix, not necessarily a directory."""
    candidates = [path] + [Path(str(path) + suffix) for suffix in (".1.bt2", ".1.bt2l", ".rev.1.bt2", ".rev.1.bt2l")]
    if any(candidate.exists() for candidate in candidates):
        log(f"OK rRNA Bowtie2 index: {path}")
    else:
        errors.append(f"Missing rRNA Bowtie2 index prefix/files: {path}")


def base_config() -> dict:
    return {
        "input_dir": "",
        "output_dir": "",
        "metadata_file": "",
        "rnaseq_count_matrix": "",
        "rnaseq_heatmap_top": "50,100",
        "comparative_group_col": "group_all",
        "comparative_ref_level": "",
        "rscript_bin": "",
        "run_profile": "full",
        "test_sample_count": 2,
        "test_read_count": 10000,
        "test_seed": 7,
        "threads": 16,
        "salmon_threads": 8,
        "rrna_index": "",
        "eggnog_data_dir": "",
        "fastp_strategy": "full",
        "fastp_phred": 20,
        "fastp_min_len": 50,
        "fastp_extra": "--detect_adapter_for_pe",
        "megahit_k_list": "21,41,61,81,101,121",
        "megahit_min_contig_len": 300,
        "quast_mode": "quast",
        "tpm_threshold": 1.0,
        "tpm_min_samples": 1,
        "card_tpm_threshold": 1.0,
        "card_cutoffs": "Strict,Loose",
        "resume": True,
        "force": False,
        "start_from": "beginning",
        "recursive": True,
        "steps": [
            "estimate_rrna", "trim_qc", "remove_rrna", "rrna_stats", "megahit",
            "collect_contigs", "quast", "transdecoder", "collect_cds", "cds_stats",
            "salmon", "salmon_matrix", "eggnog", "aggregate", "kegg_pathways",
            "comparative", "rnaseq_overview", "amr_ko", "card_rgi", "card_summary",
            "card_ko_integration", "multiqc", "summary",
        ],
    }


def write_config(prefix: Path, envs: dict[str, str], databases: dict[str, str], engine: str, runtime: Path) -> Path:
    config = base_config()
    config["envs"] = envs
    config["rrna_index"] = databases["rrna_index"]
    config["eggnog_data_dir"] = databases["eggnog_data_dir"]
    if databases.get("card_data_dir"):
        config["card_data_dir"] = databases["card_data_dir"]
    config["conda_bin"] = engine
    config["rscript_bin"] = str(Path(envs["r"]) / "bin" / "Rscript")
    config["installer_prefix"] = str(prefix)
    config["runtime_dir"] = str(runtime)
    output = prefix / "metarex_config.json"
    output.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return output


def install(args: argparse.Namespace) -> int:
    if platform.system() != "Linux":
        raise RuntimeError("The supplied lockfiles are linux-64 lockfiles; install on Linux or provide new locks")
    lock_dir = find_lock_dir()

    prefix = Path(args.prefix).expanduser().resolve()
    env_root = Path(args.env_root).expanduser().resolve()
    db_root = prefix / "databases"
    engine = find_engine(args.engine)
    log(f"Using package manager: {engine}")
    envs: dict[str, str] = {}
    environment_names = {
        "qc": "metagenomics_base",
        "meta": "metatrascriptomics_base",
        "quast": "nanopore_assembly",
        "bioinfo": "bioinfo_base",
        "eggnog": "eggnog_mapper_v2",
        "card": "card",
    }
    for key, lock_name in LOCKS.items():
        env_path = Path(envs["meta"]) if key == "r" and "meta" in envs else env_root / environment_names[key]
        envs[key] = str(env_path)
    for key, lock_name in LOCKS.items():
        if key == "r":
            continue
        engine_create(engine, Path(envs[key]), lock_dir / lock_name, prefix / "pkgs")

    # The R environment intentionally shares the meta environment.  Install
    # the R packages once, after the explicit lockfile has created that env.
    install_r_packages(engine, Path(envs["r"]), prefix / "pkgs")

    databases = {
        "rrna_index": str(Path(args.rrna_index).expanduser().resolve()) if args.rrna_index else str(db_root / "SILVA_138_2_rRNA"),
        "eggnog_data_dir": str(Path(args.eggnog_data).expanduser().resolve()) if args.eggnog_data else str(db_root / "eggnog"),
    }
    if args.card_data:
        databases["card_data_dir"] = str(Path(args.card_data).expanduser().resolve())
    if args.rrna_url:
        databases["rrna_index"] = str(download_and_extract(args.rrna_url, db_root / "SILVA"))
    if args.eggnog_url:
        databases["eggnog_data_dir"] = str(download_and_extract(args.eggnog_url, db_root / "eggNOG"))
    if args.card_url:
        databases["card_data_dir"] = str(download_and_extract(args.card_url, db_root / "CARD"))

    prefix.mkdir(parents=True, exist_ok=True)
    runtime = copy_runtime(prefix)
    config_path = write_config(prefix, envs, databases, engine, runtime)
    report = {
        "platform": platform.platform(),
        "package_manager": engine,
        "prefix": str(prefix),
        "environments": envs,
        "databases": databases,
        "runtime_dir": str(runtime),
    }
    (prefix / "installation.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    log(f"Configuration written to {config_path}")
    return check_config(config_path, engine, allow_missing_databases=True)


def load_config(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def check_config(config_path: Path, engine: str | None = None, *, allow_missing_databases: bool = False) -> int:
    config = load_config(config_path)
    engine = engine or find_engine("auto")
    errors: list[str] = []
    for key, env in config.get("envs", {}).items():
        check_path(f"environment {key}", Path(env), errors)
        for tool in TOOLS.get(key, ()):
            result = run_in_env(engine, Path(env), ["which", tool])
            if result.returncode == 0:
                log(f"OK {key}/{tool}")
            else:
                errors.append(f"Missing executable {tool} in {env}")

    r_env_text = config.get("envs", {}).get("r") or config.get("envs", {}).get("meta")
    if r_env_text:
        r_env = Path(r_env_text)
        for r_package in R_PACKAGES:
            expression = f"if (!requireNamespace('{r_package}', quietly=TRUE)) quit(status=1)"
            result = run_in_env(engine, r_env, ["Rscript", "-e", expression])
            if result.returncode == 0:
                log(f"OK R/{r_package}")
            else:
                errors.append(f"Missing R package {r_package} in {r_env}")
        meta_env_text = config.get("envs", {}).get("meta")
        if meta_env_text:
            meta_env = Path(meta_env_text)
            for import_name in PYTHON_PACKAGES:
                result = run_in_env(
                    engine,
                    meta_env,
                    ["python", "-c", f"import importlib; importlib.import_module('{import_name}')"],
                )
                if result.returncode == 0:
                    log(f"OK Python/{import_name}")
                else:
                    errors.append(f"Missing Python package {import_name} in {meta_env}")
    database_errors: list[str] = []
    check_rrna_index(Path(config.get("rrna_index", "")), database_errors)
    check_path("eggNOG data directory", Path(config.get("eggnog_data_dir", "")), database_errors)
    if allow_missing_databases:
        for warning in database_errors:
            log("WARNING " + warning)
    else:
        errors.extend(database_errors)
    report = config_path.parent / "installation_check.tsv"
    report.write_text("status\tmessage\n" + "\n".join(
        ("FAIL\t" + error for error in errors)
    ) + ("\n" if errors else "OK\tall checks passed\n"), encoding="utf-8")
    if errors:
        for error in errors:
            log("ERROR " + error)
        return 1
    log("All MetaREX installation checks passed")
    return 0


def versions(config_path: Path) -> int:
    config = load_config(config_path)
    engine = find_engine("auto")
    rows = ["environment\ttool\tversion"]
    for key, env in config.get("envs", {}).items():
        for tool in TOOLS.get(key, ()):
            result = run_in_env(engine, Path(env), [tool, "--version"])
            output = (result.stdout + result.stderr).strip().replace("\n", " ")
            rows.append(f"{key}\t{tool}\t{output or 'unavailable'}")
    output_path = config_path.parent / "software_versions.tsv"
    output_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(output_path)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Install and validate MetaREX")
    sub = parser.add_subparsers(dest="command", required=True)
    install_parser = sub.add_parser("install", help="create locked environments and write configuration")
    install_parser.add_argument("--prefix", default=str(Path.home() / "metarex"))
    install_parser.add_argument(
        "--env-root",
        default="/data/conda_envs",
        help="directory containing the script-compatible Conda environments (default: /data/conda_envs)",
    )
    install_parser.add_argument("--engine", choices=("auto", "conda", "mamba", "micromamba"), default="auto")
    install_parser.add_argument("--rrna-index", help="existing Bowtie2 SILVA index prefix")
    install_parser.add_argument("--eggnog-data", help="existing eggNOG-mapper data directory")
    install_parser.add_argument("--card-data", help="reserved for an existing CARD database path")
    install_parser.add_argument("--rrna-url", help="explicit SILVA archive URL")
    install_parser.add_argument("--eggnog-url", help="explicit eggNOG archive URL")
    install_parser.add_argument("--card-url", help="explicit CARD archive URL")
    check_parser = sub.add_parser("check", help="validate an installed configuration")
    check_parser.add_argument("--config", required=True, type=Path)
    versions_parser = sub.add_parser("versions", help="write software_versions.tsv")
    versions_parser.add_argument("--config", required=True, type=Path)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.command == "install":
        return install(args)
    if args.command == "check":
        return check_config(args.config)
    return versions(args.config)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, OSError, subprocess.CalledProcessError) as error:
        print(f"[MetaREX installer] ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
