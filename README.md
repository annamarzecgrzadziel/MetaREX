# MetaREX: RNA-seq metatranscriptomic AMR pipeline

MetaREX is a local HTML application for quantitative analysis of active
antimicrobial-resistance signals in paired-end environmental RNA-seq data.
The repository contains the application, supporting scripts, and the numbered
source scripts in [`scripts/`](scripts/).

## Source-code provenance

The files in `scripts/` are the final, frozen source-script set used for the
reported pipeline analysis. They are distributed without modification as part
of the public release. The run configuration and pipeline log should be
retained with every analysis output.

## Requirements

- Linux with Python 3.10 or newer;
- the command-line tools and R packages listed in
  [`pipeline_publication_summary.md`](docs/pipeline_publication_summary.md);
- paired-end FASTQ files with matching R1/R2 names;
- a tab-separated metadata file whose first column contains sample names;
- a Bowtie2 SILVA rRNA index and an eggNOG-mapper data directory;
- R with DESeq2, ALDEx2, apeglm, tidyverse, ggplot2 and pheatmap.

## Quick start

Start the local interface:

```bash
./start_rnaseq_amr_pipeline_app.sh
```

Then open `http://127.0.0.1:8791` and enter the input directory, output
directory, metadata file, database paths, and resource settings.

The Conda environment paths can be supplied in the JSON configuration under
`envs`. Alternatively, set `METAREX_ENV_QC`, `METAREX_ENV_META`,
`METAREX_ENV_QUAST`, `METAREX_ENV_EGGNOG`, `METAREX_ENV_CARD`, and
`METAREX_ENV_R`. The database paths can be supplied as `rrna_index` and
`eggnog_data_dir`, or through `METAREX_RRNA_INDEX` and
`METAREX_EGGNOG_DATA`.

For a reproducible command-line run, copy
[`example_pipeline_config.json`](example_pipeline_config.json), replace the
placeholder paths, and run:

```bash
python3 rnaseq_amr_pipeline_app.py run \
  --config example_pipeline_config.json
```

The application writes `pipeline_config.json`, `pipeline.log`, intermediate
results, and a final `index.html` into the configured output directory.

The default paths shown in the graphical interface are relative placeholders;
all input, output, environment and database locations must be adapted to the
local installation.

## Expected metadata

The metadata file must contain a sample identifier column matching the FASTQ
sample names. Comparative analyses use the `group_all` column by default; this
can be changed with `comparative_group_col`.

## Main workflow

1. FASTQ discovery and optional test-run downsampling;
2. fastp quality control and MultiQC;
3. Bowtie2 rRNA estimation and removal;
4. independent per-sample MEGAHIT assembly;
5. QUAST/MetaQUAST assembly QC;
6. TransDecoder CDS and peptide prediction;
7. Salmon quantification and expression matrices;
8. eggNOG-mapper annotation and KO/EC/PFAM aggregation;
9. KO-to-KEGG Pathway aggregation;
10. ALDEx2 and DESeq2 comparative analyses;
11. VST/PCA/heatmap exploratory analysis;
12. KO-based AMR, CARD/RGI screening, and CARD–KO integration.

Cross-sample statistical analyses use common functional features obtained after
eggNOG-based aggregation. Independent sample assemblies are not treated as a
single shared transcriptome.

## Reproducibility records

The publication-release documentation records the software versions,
parameters, databases, checksums, and methodological decisions:

- [`pipeline_publication_summary.md`](docs/pipeline_publication_summary.md)
- [`reproducibility_and_methods.md`](docs/reproducibility_and_methods.md)
- [`software_versions.md`](provenance/software_versions.md)
- [`environment_report.md`](docs/environment_report.md)
- [`SHA256SUMS`](checksums/public_release_scripts_SHA256SUMS)
