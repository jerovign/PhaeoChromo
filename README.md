# PhaeoChromo
### ***Emergence and evolution of chromatin landscapes across brown algae (Phaeophyceae)***

Jérômine Vigneau<sup>1,\*</sup>, Jaruwatana Sodai Lotharukpong<sup>1,\*</sup>, Pengfei Liu<sup>1</sup>, Rémy Luthringer<sup>1</sup>, Bérangère Lombard<sup>2</sup>, Damarys Loew<sup>2</sup>, Fabian Haas<sup>1</sup>, Michael Borg<sup>1,†</sup>, Susana M. Coelho<sup>1,†</sup> 

<sup>1</sup>Max Planck Institute for Biology, Tübingen, Germany

<sup>2</sup>Institut Curie, PSL Research University, Paris, France

† Co-corresponding authors

## Overview

This repository contains the analysis code and figure generation scripts associated with the study:

> **Vigneau et al.** *Emergence and evolution of chromatin landscapes across brown algae*
(manuscript in preparation)

We investigate how chromatin landscapes emerged and diversified across multiple brown algal lineages. By combining genome-wide histone modification profiles, transcriptomics, and evolutionary analyses, we:
- identify **conserved chromatin signatures** linked to gene activity,
- uncover **lineage- and sex-specific chromatin dynamics**,
- and explore their role in **sex chromosome evolution**.

## Repository orgnazation

* `01_Make_tables_from_geneID_to_chromatin-signatures/` — Notebooks describing the processed data tables (containing list of genes, their expression, chromatin signature assignment and other details).

* `02_Make_chromatin_states_with_hiHMM_software/` — Notebooks explaining the workflow for making input files for hiHMM and running the software.

* `Scripts/` — Other R, Python, bash scripts used for specific tasks - written on the subfolder.

## Requirements

Analyses were performed in R (≥ 4.3) with:
- [`tidyverse v2.0.0`](https://www.tidyverse.org/) (incl. `ggplot2 v3.5.2`, `dplyr v1.1.4`, `tidyr v1.3.1`, `tibble v3.3.0`, `stringr v1.5.1`, `readxl v1.4.3`)
- [`rtracklayer v1.64.0`](https://www.bioconductor.org/packages/release/bioc/html/rtracklayer.html)
- [`GenomicRanges v1.56.2`](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)
- [`ggalluvial`](https://corybrunson.github.io/ggalluvial/)
- [`ggrastr v1.0.2`](https://github.com/VPetukhov/ggrastr)
- [`rsample v1.2.1`](https://rsample.tidymodels.org/) from [tidymodels](https://www.tidymodels.org/)

Preprocessing of sequencing data used:
- Nextflow (≥ 23) nf-core pipelines:
  - [`nf-core/chipseq v2.0.0`](https://nf-co.re/chipseq/2.0.0)
  - [`nf-core/rnaseq v3.12.0`](https://nf-co.re/rnaseq/3.12.0)
- Standard bioinformatics tools (gffread...)

## Citation

If you use this repository, please cite:

> **Vigneau J., Lotharukpong J.S., Liu P., Luthringer R., Lombard B., Loew D., Haas F., Borg M.†, Coelho S.M.†**
> 
> *Emergence and evolution of chromatin landscapes across brown algae.* (in preparation)

## Contact
Maintainer: Jérômine Vigneau
