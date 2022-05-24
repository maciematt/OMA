# OMA - Omics Meta-Analysis

## Installation

Install from the github repo via

```
> devtools::install_github("maciematt/OMA")
```

or from a local copy via

```
> install.packages("location/with/OMA", repos = NULL, type = "source")
```

## Usage

Applicatin of OMA occurs in three stages: 1. data preparation and checking, 2. differential gene expression analysis in each dataset, and 3. meta-analysis of the DGE results from 2.

1. Data prep for each dataset:

```
library(OMA)

# 0. Important sample filters!
sample_filters <- rlang::quo(
  grepl("skin", tissue) & (sampling_time %in% c("No Info", "W0", "baseline", "pre-treatment")) & (disease_state == "atopic dermatitis")
)

## 1. Preparing each dataset
ad_for_dge <- lapply(1:length(ad_datasets), function (x) {
  prepare_for_dge(
    dataset = ad_datasets[[x]], ## each dataset is a list where the first member is the metadata, and the second member is the tx data
    dataset_name = names(ad_datasets)[x], ## the name of the dataset
    contrast = list(variable = "sample_pathology", active = "lesional", reference = "non-lesional"),  ## this is the main contrast - "variable" is the name of the column from the meta-data file tot use for the main contrast; active is the active level of that variable; reference is the reference level
    blocks = "subject_id", scrutinize_blocks = TRUE, ## additional OPTIONAL setup - a blocking variable
    covariates = "gender", scrutinize_covariates = TRUE, ## additional OPTIONAL setup - covariates in the model
    sample_filters = sample_filters ## additional OPTIONAL setup (not recommended) - if your data doesn't come approprietly pre-filtered, some filtering can be specified here
  )
})
```

2. DGE analysis:

```
ad_dge <- lapply(ad_for_dge, run_dge)
```

3. Meta-analysis:

```
ad_ma <- run_ma(ad_dge, es_var = "coeff", parallel = TRUE)
```

That's it. The resulting list can be accessed via `ma$ma`.

