# MetaREX: RNA-seq AMR Metatranscriptomics Pipeline

Local, interactive, browser-based application for running an RNA-seq / metatranscriptomics workflow with functional annotation and antimicrobial resistance (AMR) profiling.

The application is designed to be started on a local workstation or server and controlled through a web browser. It is not intended to be deployed as a public web service. All analyses are executed locally using command-line tools, Conda environments and reference databases installed on the host machine.

---

## Overview

This project provides a lightweight web interface for a multi-step RNA-seq / metatranscriptomics pipeline. The interface helps configure input/output paths, computational parameters, metadata, selected analysis steps and test runs, while the actual analyses are executed locally.

The workflow covers:

- paired-end FASTQ discovery,
- optional rRNA content estimation,
- read trimming and QC with `fastp` and `MultiQC`,
- rRNA removal with `Bowtie2`,
- de novo transcript assembly with `MEGAHIT`,
- assembly quality assessment with `QUAST` or `MetaQUAST`,
- coding sequence prediction with `TransDecoder`,
- transcript quantification with `Salmon`,
- functional annotation with `eggNOG-mapper`,
- KO / EC / PFAM / KEGG aggregation,
- differential analysis with `ALDEx2` and `DESeq2`,
- exploratory RNA-seq plots,
- AMR-oriented KO summaries,
- CARD/RGI-based AMR profiling,
- CARD × KO integration,
- final MultiQC and HTML summary reports.

---

## Important note

This is a local analysis application.

- FASTQ files are not uploaded anywhere.
- Reference databases remain local.
- Intermediate and final results are written to the selected local output directory.
- The default web interface binds to `127.0.0.1`, so it is accessible only from the same machine unless the host is changed manually.

---

## Repository structure

Recommended repository layout:

```text
MetaREX/
├── README.md
├── .gitignore
├── rnaseq_amr_pipeline_app.py
├── start_rnaseq_amr_pipeline_app.sh
├── scripts/
│   ├── 00_estimate_rRNA_content.sh
│   ├── 01_fastp_trim_QC.sh
│   ├── 02_remove_rRNA_bowtie2.sh
│   ├── 02a_remove_rRNA_stats.sh
│   ├── 03a_megahit.sh
│   ├── 03b_contigs_all.sh
│   ├── 03c_quast_metaquast.sh
│   ├── 04_transdecoder.sh
│   ├── 04a_transdecoder_all_cds.sh
│   ├── 05_compare_cds_stats.py
│   ├── 06_quant_salmon.sh
│   ├── 07_collect_salmon_tpm.py
│   ├── 08_eggnog_mapper2.sh
│   ├── 09a_aggregate_tpm_by_eggnog.py
│   ├── 09b_aggregate_counts_by_eggnog.py
│   ├── 10a_DEG_KO_ALDEx2.R
│   ├── 10a_DEG_KEGG_Pathway_ALDEx2.R
│   ├── 10b_DEG_KO_DeSeq2.R
│   ├── 10b_DEG_KEGG_Pathway_DeSeq2.R
│   ├── 11_KO_to_KEGG_Pathway.R
│   ├── 12_KEGG_pathway_categories.R
│   ├── 13_RNAseq_one_group_analysis.R
│   ├── 14_AMR_KO.R
│   ├── 15_CARD_RGI_protein.sh
│   ├── 15a_CARD_analysis.R
│   └── 16_CARD_KO_integration.R
│  
├── docs/
│   └── workflow_overview.md
├── results/
│   └── .gitkeep
└── logs/
    └── .gitkeep
```

Large input files, databases and generated results should not be committed to Git.

---

## Requirements

### Operating system

The pipeline is intended for Linux-based systems, for example:

- Ubuntu,
- Debian,
- CentOS / Rocky / Oracle Linux,
- a Linux server accessed through SSH,
- WSL2 on Windows, if all external tools and databases are available.

### Python

The web application itself uses only Python standard library modules.

Recommended:

```bash
python3 --version
```

Python 3.10 or newer is recommended.

### External command-line tools

The pipeline expects the following tools to be installed in local Conda environments or available in `PATH`:

- `fastp`
- `FastQC`
- `MultiQC`
- `Bowtie2`
- `MEGAHIT`
- `QUAST` / `MetaQUAST`
- `TransDecoder`
- `Salmon`
- `eggNOG-mapper`
- `seqtk`
- `RGI` from CARD
- `Rscript`

### R packages

Depending on selected steps, the following R packages may be required:

- `tidyverse`
- `DESeq2`
- `ALDEx2`
- `ggplot2`
- `pheatmap`
- `apeglm`
- `matrixStats`
- `janitor`
- `FactoMineR`
- `factoextra`
- `KEGGREST`

### Databases and indexes

The following local resources are expected:

- Bowtie2 rRNA index, for example SILVA rRNA,
- eggNOG-mapper data directory,
- CARD/RGI database,
- any additional reference resources required by selected tools.

Example local paths used during development:

```text
/user_bases_localization/SILVA_138_2_rRNA
/user_bases_localization/eggnog
/user_conda_envs_localization/metagenomics_base
/user_conda_envs_localization/metatrascriptomics_base
/user_conda_envs_localization/eggnog_mapper_v2
/user_conda_envs_localization/card
```

These paths should be adapted to the local machine before running the pipeline.

---

## Installation

Clone the repository:

```bash
git clone https://github.com/YOUR-USERNAME/MetaREX.git
cd MetaREX
```

Make the launcher executable:

```bash
chmod +x start_rnaseq_amr_pipeline_app.sh
```

Check that Python is available:

```bash
python3 --version
```

Check that the required Conda environments and databases are available on the local machine.

---

## Local configuration

Before the first run, review paths in:

```text
rnaseq_amr_pipeline_app.py
```

In particular, check the following sections:

```python
ENVS = {
    "qc": "/user_conda_envs_localization/metagenomics_base",
    "meta": "/user_conda_envs_localization/metatrascriptomics_base",
    "quast": "/user_conda_envs_localization/nanopore_assembly",
    "eggnog": "/user_conda_envs_localization/eggnog_mapper_v2",
    "card": "/user_conda_envs_localization/card",
    "r": "/user_conda_envs_localization/metatrascriptomics_base",
}

DBS = {
    "rrna_bowtie2": "/user_bases_localization/SILVA_138_2_rRNA",
    "eggnog_data": "/user_bases_localization/eggnog",
}
```

If your scripts are stored inside the repository under `scripts/`, update the `SOURCE_SCRIPTS` paths accordingly, for example:

```python
SOURCE_SCRIPTS = {
    "fastp": "scripts/01_fastp_trim_QC.sh",
    "remove_rrna": "scripts/02_remove_rRNA_bowtie2.sh",
    "megahit": "scripts/03a_megahit.sh",
    # ...
}
```

Alternatively, keep absolute paths if the pipeline is used only on a dedicated local server.

---

## Running the browser application

Start the local web app:

```bash
bash start_rnaseq_amr_pipeline_app.sh
```

By default, the app starts at:

```text
http://127.0.0.1:8791
```

Open this address in a web browser on the same machine.

You can also start the app directly:

```bash
python3 rnaseq_amr_pipeline_app.py serve --host 127.0.0.1 --port 8791
```

To change the host or port:

```bash
HOST=127.0.0.1 PORT=8792 bash start_rnaseq_amr_pipeline_app.sh
```

If the app is running on a remote server, use SSH port forwarding rather than exposing the app publicly:

```bash
ssh -L 8791:127.0.0.1:8791 user@server
```

Then open locally:

```text
http://127.0.0.1:8791
```



## Input FASTQ files

The input directory should contain paired-end FASTQ files, preferably gzipped.

Supported extensions include:

```text
.fastq
.fastq.gz
.fq
.fq.gz
```

The application automatically searches for paired-end files using common R1/R2 naming patterns, for example:

```text
sample_01_R1.fastq.gz
sample_01_R2.fastq.gz
```

or Illumina-style names such as:

```text
sample_01_S1_L001_R1_001.fastq.gz
sample_01_S1_L001_R2_001.fastq.gz
```

The app can search recursively inside the selected input directory. Derived FASTQ files such as cleaned, rRNA-free or test-run files are ignored where possible to avoid accidental reanalysis of intermediate data.

---

## Metadata file

A metadata file is optional but recommended for comparative analyses and plotting.

Supported formats:

- TSV
- CSV

The metadata file should contain at least one of the following columns:

```text
sample
```

or:

```text
sampleID
```

Recommended minimal TSV format:

```text
sampleID	group
sample_01	control
sample_02	control
sample_03	treated
sample_04	treated
```

Optional columns such as `batch`, `condition`, `site`, `timepoint` or other experimental variables can be added if used by downstream analyses.

If no metadata file is provided, the application creates a default metadata table in the output directory and assigns all samples to a single group called `all`.

---

## Main configuration fields in the web interface

### Paths

| Field | Description |
|---|---|
| `Input FASTQ directory` | Directory containing paired-end FASTQ files. |
| `Output directory` | Directory where all pipeline results will be written. |
| `Metadata TSV/CSV` | Optional metadata file. |

### RNA-seq and comparisons

| Field | Description |
|---|---|
| `RNA-seq count matrix` | Optional external count matrix for overview plots. |
| `Heatmap top genes` | Comma-separated numbers of top features for heatmaps, for example `50,100`. |
| `Grouping column` | Metadata column used for group-based comparisons. |
| `Reference group` | Reference level for DESeq2-based comparisons. |
| `Rscript binary` | Path to the `Rscript` executable. |

### Run profile

| Option | Description |
|---|---|
| `Full run` | Uses all discovered reads. |
| `Test run` | Downsamples selected samples to a smaller number of reads for testing. |

### Resources and databases

| Field | Description |
|---|---|
| `Threads` | General number of CPU threads. |
| `Salmon threads` | Threads used by Salmon quantification. |
| `Bowtie2 rRNA index` | Prefix of the Bowtie2 rRNA index. |
| `eggNOG data directory` | Local eggNOG-mapper data directory. |

### Parameters

| Field | Description |
|---|---|
| `fastp strategy` | Predefined or custom trimming strategy. |
| `QUAST mode` | `QUAST` or `MetaQUAST`. |
| `fastp Phred` | Minimum quality threshold for selected fastp modes. |
| `fastp min length` | Minimum read length after trimming. |
| `Custom fastp options` | Custom fastp parameters when custom mode is selected. |
| `MEGAHIT k-list` | Comma-separated k-mer list. |
| `MEGAHIT min contig length` | Minimum contig length. |
| `TPM threshold` | Expression threshold for TPM filtering. |
| `Min samples for TPM filter` | Minimum number of samples meeting TPM threshold. |
| `CARD TPM threshold` | TPM threshold for CARD/RGI analysis on expressed contigs. |
| `CARD cutoffs` | CARD/RGI cutoffs to retain, for example `Strict,Loose`. |

### Execution scope

| Field | Description |
|---|---|
| `Start from step` | Allows restarting the pipeline from a selected stage. |
| `Run mode` | `resume` skips completed steps; `force` reruns selected steps. |
| `Steps` | Individual workflow steps to include or exclude. |

---

## Workflow steps

| Step ID | Description |
|---|---|
| `estimate_rrna` | Estimate rRNA content using Bowtie2. |
| `trim_qc` | Trim reads and generate fastp/MultiQC reports. |
| `remove_rrna` | Remove rRNA reads using Bowtie2. |
| `rrna_stats` | Summarize read retention after rRNA removal. |
| `megahit` | Assemble rRNA-free reads with MEGAHIT. |
| `collect_contigs` | Collect final contigs from sample-level assemblies. |
| `quast` | Generate QUAST or MetaQUAST assembly reports. |
| `transdecoder` | Predict coding sequences and proteins with TransDecoder. |
| `collect_cds` | Collect CDS and protein FASTA files. |
| `cds_stats` | Summarize CDS statistics. |
| `salmon` | Quantify expression with Salmon. |
| `salmon_matrix` | Generate TPM and count matrices. |
| `eggnog` | Annotate predicted proteins with eggNOG-mapper. |
| `aggregate` | Aggregate expression/counts by KO, EC and PFAM. |
| `kegg_pathways` | Convert KO-level data to KEGG pathway-level matrices. |
| `comparative` | Run group comparisons using ALDEx2 and DESeq2. |
| `rnaseq_overview` | Generate exploratory RNA-seq plots. |
| `amr_ko` | Summarize AMR-related KO expression. |
| `card_rgi` | Run CARD/RGI on expressed contigs/proteins. |
| `card_summary` | Summarize CARD/RGI results. |
| `card_ko_integration` | Integrate CARD and KO-derived AMR signals. |
| `multiqc` | Generate final MultiQC report. |
| `summary` | Generate final HTML index for the run. |

---

## Output structure

A typical output directory contains:

```text
output_run/
├── pipeline_config.json
├── pipeline.log
├── 00_samples.tsv
├── 00_metadata.tsv
├── 00_rRNA_content/
├── 01_fastp_trim_QC/
├── 02_rRNA_removed/
├── 02a_rRNA_stats/
├── 03_assembly_megahit/
├── 03b_contigs_all/
├── 03c_quast_metaquast/
├── 04_transdecoder/
├── 04a_all_cds/
├── 04b_all_pep/
├── 05_cds_stats/
├── 06_quant_salmon/
├── 07_salmon_matrices/
├── 08_eggnog/
├── 09a_aggregate_tpm_by_eggnog/
├── 09b_aggregate_counts_by_eggnog/
├── 10_comparative_analysis/
├── 11_KO_to_KEGG_Pathway/
├── 13_RNAseq_one_group/
├── 14_AMR_KO/
├── 15_CARD/
├── 15a_CARD_analysis/
├── 16_integracja_CARD_KO/
├── 99_multiqc_all/
└── index.html
```

Important output files include:

| File | Description |
|---|---|
| `pipeline.log` | Full execution log. |
| `pipeline_config.json` | Configuration used for the run. |
| `00_samples.tsv` | Detected paired-end samples. |
| `00_metadata.tsv` | Metadata used by the pipeline. |
| `07_salmon_matrices/salmon_TPM_all.tsv` | Transcript-level TPM matrix. |
| `09b_aggregate_counts_by_eggnog/09b_COUNTS_KO_matrix.tsv` | KO count matrix. |
| `11_KO_to_KEGG_Pathway/11_COUNTS_KEGG_Pathway_matrix.tsv` | KEGG pathway count matrix. |
| `10_comparative_analysis/` | ALDEx2 and DESeq2 comparison results. |
| `14_AMR_KO/` | AMR-related KO summaries. |
| `15_CARD/` | Per-sample CARD/RGI outputs. |
| `16_integracja_CARD_KO/` | Integrated CARD × KO outputs. |
| `99_multiqc_all/multiqc_all.html` | Final MultiQC report. |
| `index.html` | Final HTML run summary. |

---

## Test run mode

Test mode is useful for checking whether paths, tools and databases are configured correctly before running a full analysis.

In test mode, the application:

- selects a subset of samples,
- randomly downsamples paired FASTQ files,
- writes test FASTQ files to the output directory,
- runs selected steps on the reduced dataset.

Recommended first test:

```text
Test samples: 1-2
Reads per sample: 10000
Random seed: 7
```

Once the test run finishes successfully, switch to `Full run`.

---

## Resume and rerun behavior

The application records step status in the output directory and can skip steps that have already produced expected output files.

Use:

```text
Run mode: resume
```

for normal continuation of an interrupted or partially completed run.

Use:

```text
Run mode: force
```

when you want to rerun selected steps even if outputs already exist.

The `Start from step` option allows restarting the workflow from a selected stage, for example from `salmon`, `eggnog` or `comparative`.

---

## Recommended first run

For a new installation, start with a small test run:

1. Start the app:

   ```bash
   bash start_rnaseq_amr_pipeline_app.sh
   ```

2. Open:

   ```text
   http://127.0.0.1:8791
   ```

3. Select:

   ```text
   Run profile: Test run
   Test samples: 1
   Reads per sample: 10000
   ```

4. Provide:

   - input FASTQ directory,
   - output directory,
   - Bowtie2 rRNA index,
   - eggNOG data directory,
   - metadata file, if available.

5. Run only early steps first:

   ```text
   trim_qc
   remove_rrna
   rrna_stats
   megahit
   ```

6. If the test succeeds, run the full selected workflow.

---

## Data and Git policy

Do not commit large sequencing files or generated outputs.

Recommended `.gitignore` entries:

```gitignore
*.fastq
*.fastq.gz
*.fq
*.fq.gz
*.bam
*.sam
*.cram
*.bai
*.crai
*.fa
*.fasta
*.fna
*.gtf
*.gff
*.gff3
*.bt2
*.bt2l
*.idx
*.dmnd
results/*
logs/*
runs/*
work/*
tmp/*
*.log
.env
config/local*
```

Keep only:

- source code,
- scripts,
- small example metadata/config files,
- documentation,
- small test fixtures if legally and ethically shareable.

---

## Troubleshooting

### The browser does not open automatically

The launcher starts the local server but may not open the browser automatically. Open manually:

```text
http://127.0.0.1:8791
```

### No FASTQ pairs detected

Check that input files use recognizable R1/R2 naming, for example:

```text
sample_R1.fastq.gz
sample_R2.fastq.gz
```

Also check whether recursive search is enabled.

### Bowtie2 rRNA index not found

Make sure the value in `Bowtie2 rRNA index` points to the index prefix, not to a single `.bt2` file.

For example, if the directory contains:

```text
SILVA_138_2_rRNA.1.bt2
SILVA_138_2_rRNA.2.bt2
SILVA_138_2_rRNA.3.bt2
SILVA_138_2_rRNA.4.bt2
```

then the correct prefix is:

```text
/path/to/SILVA_138_2_rRNA
```

### eggNOG-mapper fails

Check:

- Conda environment path,
- `emapper.py` availability,
- eggNOG database directory,
- available disk space,
- number of CPU threads.

### R-based steps fail

Check that `Rscript` points to the correct R installation and that all required R packages are installed in that environment.

### CARD/RGI fails

Check:

- CARD Conda environment,
- `rgi` availability,
- CARD database installation,
- `seqtk` availability if sequence extraction is required.

### A step is skipped unexpectedly

Check:

- selected checkboxes in the web interface,
- `Start from step`,
- `Run mode`,
- existing output files in the run directory,
- `.pipeline_state` inside the output directory.

Use `force` mode if a step needs to be rerun.

---

## Security notes

By default, the application binds to:

```text
127.0.0.1
```

This is appropriate for local use.

Do not expose the application directly to the public internet. If remote access is needed, use SSH port forwarding or a secured internal network.

---

## Limitations

- The pipeline assumes paired-end RNA-seq/metatranscriptomic data.
- The workflow depends on local Conda environments and local databases.
- Absolute paths may need to be adapted before running on another machine.
- Some analysis steps require sufficient RAM and disk space, especially assembly, annotation and CARD/RGI steps.
- The biological interpretation of AMR-related KO and CARD/RGI results should be validated in the context of the dataset and experimental design.

---

## Citation

If you use this pipeline in a report, manuscript or presentation, please cite the underlying tools used in the selected workflow, including `fastp`, `Bowtie2`, `MEGAHIT`, `QUAST`, `TransDecoder`, `Salmon`, `eggNOG-mapper`, `ALDEx2`, `DESeq2`, `CARD` and `RGI`, as appropriate.

---

## Author

Maintainer: Anna Marzec-Grzadziel, Institute of Soil Science and Plant Cultivation, agrzadziel@iung.pulawy.pl

Project type: local RNA-seq / metatranscriptomics / AMR analysis application.
