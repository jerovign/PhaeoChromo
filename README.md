# PhaeoChromo
### ***Emergence and evolution of chromatin landscapes across brown algae (Phaeophyceae)***

Jérômine Vigneau<sup>1,*</sup>, Jaruwatana Sodai Lotharukpong<sup>1,*</sup>, Pengfei Liu<sup>1</sup>, Rémy Luthringer<sup>1</sup>, Bérangère Lombard<sup>2</sup>, Damarys Loew<sup>2</sup>, Fabian Haas<sup>1</sup>, Michael Borg<sup>1,†</sup>, Susana M. Coelho<sup>1,†</sup> 

<sup>1</sup>Max Planck Institute for Biology, Tübingen, Germany

<sup>2</sup>Institut Curie, PSL Research University, Paris, France

† Co-corresponding authors

## Overview

This repository contains the analysis code and figure generation scripts associated with the study:

**Vigneau J., Lotharukpong J.S., Liu P., Luthringer R., Lombard B., Loew D., Haas F., Borg M.†, Coelho S.M.†**

*Emergence and evolution of chromatin landscapes across brown algae*
(manuscript in preparation)

We investigate how chromatin landscapes emerged and diversified across multiple brown algal lineages. By combining genome-wide histone modification profiles, transcriptomics, and evolutionary analyses, we:
- identify **conserved chromatin signatures** linked to gene activity,
- uncover **lineage- and sex-specific chromatin dynamics**,
- and explore their role in **sex chromosome evolution**.

## Repository contents

* `Scripts/` — R, Python, bash scripts for data processing, chromatin signature assignment, statistical analyses and plots.

* `Data/` — Processed data tables (containing list of genes, their expression, signatures and other details).

* `Notebooks/` — Notebooks explaining the workflows for data processing.

## Requirements

Analyses were performed in R (≥ 4.3) with:
- `tidyverse`
- `ggplot2`, `ggalluvial`
Preprocessing of sequencing data used:
- Nextflow (≥ 23)
- Standard bioinformatics tools (gffread...)

## Citation

If you use this repository, please cite:
**Vigneau J., Lotharukpong J.S., Liu P., Luthringer R., Lombard B., Loew D., Haas F., Borg M.†, Coelho S.M.†**

*Emergence and evolution of chromatin landscapes across brown algae.* (in preparation)

## Contact
Maintainer: Jérômine Vigneau
