## Versions recorded for the canonical analysis

| Component | Version | Evidence |
|---|---:|---|
| Python | 3.10.12 | publication provenance |
| R | 4.5.3 | publication provenance; requires final run verification |
| fastp | 1.0.1 | publication provenance and local audit |
| MultiQC | 1.34 | publication provenance |
| Bowtie2 | 2.5.4 | publication provenance and local audit |
| MEGAHIT | 1.2.9 | publication provenance and local audit |
| QUAST/MetaQUAST | 5.3.0 | publication provenance and local audit |
| TransDecoder | 5.7.1 | publication provenance and local audit |
| Salmon | 1.10.3 | publication provenance and local audit |
| eggNOG-mapper | 2.1.13 | publication provenance and local audit |
| DIAMOND | 2.0.15 | eggNOG environment provenance |
| MMseqs2 | 16.747c6 | eggNOG environment provenance |
| RGI | 3.2.1 | publication provenance; audit `rgi -sv` |
| CARD database | 1.1.9 | audit `rgi -dv` |

## R packages recorded for the analysis

DESeq2 1.48.2, ALDEx2 1.40.0, apeglm 1.30.0, tidyverse 2.0.0,
ggplot2 4.0.0, dplyr 1.1.4, readr 2.1.5, stringr 1.5.2,
matrixStats 1.5.0, pheatmap 1.0.13, janitor 2.2.1, FactoMineR 2.13,
factoextra 1.0.7, Biobase 2.68.0, S4Vectors 0.46.0 and
SummarizedExperiment 1.38.1.

## Reference resources

- SILVA rRNA index: SILVA 138.2, Bowtie2 index prefix;
- eggNOG data: expected eggNOG database version 5.0.2;
- CARD/RGI data: RGI database version 1.1.9.

The eggNOG installation currently emits a warning and reports its installed
database version as unknown. This is an outstanding verification item for the
exact publication environment and should be resolved before a final public
release claim.

## Lockfiles and checksums

- Conda lockfiles for the available audit environments are in
  [`../environment_lockfiles/`](../environment_lockfiles/);
- reference database checksums are in
  [`../checksums/reference_databases_SHA256SUMS`](../checksums/reference_databases_SHA256SUMS);
- source and supporting-script checksums are in
  [`../checksums/public_release_scripts_SHA256SUMS`](../checksums/public_release_scripts_SHA256SUMS);
- the metadata checksum is in
  [`../checksums/example_input_SHA256SUMS`](../checksums/example_input_SHA256SUMS).

