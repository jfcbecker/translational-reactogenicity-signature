# translational-reactogenicity-signature

Code accompanying the publication:

> **A translational transcriptomic signature of vaccine reactogenicity for the evaluation of novel formulations**
> Jérémie Becker, Maroussia Roelens, Kendra Reynaud, Laurent Beloeil (BIOASTER, Lyon, France)
> *eLife* (reviewed preprint), 2026. DOI: [10.7554/eLife.109928.1](https://doi.org/10.7554/eLife.109928.1)

## Overview

This repository contains the R analysis pipeline used to derive and validate a cross-species transcriptomic signature of vaccine reactogenicity. A penalized ordinal regression model is trained on mouse muscle transcriptomes from seven vaccines and immunostimulants, then transferred to mouse blood and human blood. The resulting signature highlights inflammatory programs (notably IL-6 / JAK / STAT3) and correctly ranks licensed human vaccines by their clinical reactogenicity, with the MF59-adjuvanted influenza vaccine identified as the most reactogenic formulation.

## Data

Transcriptomic datasets are publicly available from the NCBI Gene Expression Omnibus, originating from the BioVacSafe consortium:

- [GSE120661](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120661)
- [GSE124533](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124533)

Place the downloaded data under the path defined by `dataPath` in [sibylle_mouse_human_data_preprocessing_parameters.yaml](sibylle_mouse_human_data_preprocessing_parameters.yaml).

## Repository contents

| File | Purpose |
| --- | --- |
| [sibylle_parameters.R](sibylle_parameters.R) | Global parameters (seeds, thresholds, palettes, reactogenicity classes) |
| [sibylle_mouse_human_data_preprocessing_parameters.yaml](sibylle_mouse_human_data_preprocessing_parameters.yaml) | Paths, tissues, treatment colors |
| [sibylle_mouse_human_data_preprocessing_launcher.R](sibylle_mouse_human_data_preprocessing_launcher.R) | Entry point that runs all preprocessing steps |
| [sibylle_data_loader.R](sibylle_data_loader.R) | Loads counts and metadata |
| [sibylle_mouse_data_preprocessing.R](sibylle_mouse_data_preprocessing.R) | Mouse counts QC, filtering, normalization |
| [sibylle_human_data_preprocessing.R](sibylle_human_data_preprocessing.R) | Human counts QC, filtering, normalization |
| [sibylle_mouse_human_orthology_preamble.R](sibylle_mouse_human_orthology_preamble.R) / [sibylle_mouse_human_orthology.R](sibylle_mouse_human_orthology.R) | Mouse-human ortholog mapping |
| [sibylle_D0_corrections_PCA.R](sibylle_D0_corrections_PCA.R) | Day-0 baseline correction and PCA diagnostics |
| [sibylle_mouse_redo_differential_analysis.R](sibylle_mouse_redo_differential_analysis.R) | Differential expression analysis (mouse) |
| [sibylle_human_reactogenicity_readouts_analysis.R](sibylle_human_reactogenicity_readouts_analysis.R) | Human reactogenicity readouts and signature evaluation |
| [sibylle_signature_varImportance.R](sibylle_signature_varImportance.R) | Variable importance for the signature genes |
| [sibylle_mouse_human_additional_illustrations_for_paper.R](sibylle_mouse_human_additional_illustrations_for_paper.R) | Supplementary figure generation |
| [sibylle_supplementary_tables_for_paper.R](sibylle_supplementary_tables_for_paper.R) | Supplementary tables generation |
| [sibylle_functions.R](sibylle_functions.R) | Shared utility functions |

## Requirements

- R (>= 4.2 recommended)
- CRAN packages: `caret`, `cluster`, `dplyr`, `ggplot2`, `ggh4x`, `ggnewscale`, `ggpubr`, `ggvenn`, `glmnetcr`, `matrixStats`, `mclust`, `mltools`, `openxlsx`, `patchwork`, `RColorBrewer`, `reshape2`, `rsample`, `stringr`, `this.path`, `tidyr`, `tidyverse`, `yaml`
- Bioconductor packages: `clusterProfiler`, `edgeR`, `GEOquery`, `limma`, `msigdbr`, `preprocessCore`, `sva`

Install with:

```r
install.packages(c(
  "caret", "cluster", "dplyr", "ggplot2", "ggh4x", "ggnewscale",
  "ggpubr", "ggvenn", "glmnetcr", "matrixStats", "mclust", "mltools",
  "openxlsx", "patchwork", "RColorBrewer", "reshape2", "rsample",
  "stringr", "this.path", "tidyr", "tidyverse", "yaml"
))

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "clusterProfiler", "edgeR", "GEOquery", "limma",
  "msigdbr", "preprocessCore", "sva"
))
```

## Usage

1. Edit `projectPath`, `dataPath`, and `outputPath` in [sibylle_mouse_human_data_preprocessing_parameters.yaml](sibylle_mouse_human_data_preprocessing_parameters.yaml) to match your environment.
2. Run preprocessing and orthology mapping:
   ```r
   source("sibylle_mouse_human_data_preprocessing_launcher.R")
   ```
3. Run downstream analyses (differential expression, signature training, evaluation, figures) by sourcing the corresponding scripts listed above.

The pipeline parallelizes across `nbcores` cores (default 80 in [sibylle_parameters.R](sibylle_parameters.R:18)); reduce this value to suit your machine.

## Citation

If you use this code, please cite the publication above. A machine-readable citation is provided in [CITATION.cff](CITATION.cff).

## License

Released under the [MIT License](LICENSE).
