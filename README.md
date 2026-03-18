# MetaScreener: a robust dual-mode framework for directional prioritization of actionable signatures through multi-dataset and multi-approach integration

## Overview

MetaScreener is an R package that provides a comprehensive framework for screening directional gene signatures through multi-level meta-analysis. It integrates multiple enrichment methods with various statistical approaches to identify robust activation and inhibition patterns across multiple datasets.

The package offers two complementary analysis modes:

1. **DiffMetaScreener**: Based on differential expression analysis between groups
2. **CorMetaScreener**: Based on correlation analysis with continuous variables

## Installation

```r
# Clear environment and free memory
rm(list=ls())
gc()

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install required dependencies
depens = c("singscore","UCell","decoupleR","AUCell","metap",
           "fgsea","plyr","reshape2","GSEABase","GSVA","doBy",
           "clusterProfiler","limma","dplyr","tidyverse","devtools","stringr",
           "GSA","coin","viper","caret","wdm","correlation")

# Install IOBR from GitHub
devtools::install_github("IOBR/IOBR")

# Install Bioconductor packages
for(i in 1:length(depens)){
  depen = depens[i]
  if (!requireNamespace(depen, quietly = TRUE))
    BiocManager::install(depen,update = FALSE)
}

# Install MetaScreener
devtools::install_github("FangZY-Lab/MetaScreener")

# Load the package
library(MetaScreener)
```

## Data Preparation

### Required Data Formats

#### For DiffMetaScreener:
- **Expression data**: Data frames named with pattern `diff{GSEID}` containing gene expression matrices (genes as rows, samples as columns)
- **Group information**: Data frames named with pattern `diff{GSEID}_G` containing sample grouping information
- **Comparison file**: A data frame specifying group comparisons and their biological status

#### For CorMetaScreener:
- **Expression data**: Data frames named with pattern `cor{GSEID}` containing gene expression matrices
- **Value data**: Data frames named with pattern `cor{GSEID}_V` containing continuous variables for correlation analysis

### Example Data Loading

```r
# Set working directory to your data location
setwd("path/to/your/data")
# The sample data download site is: https://github.com/FangZY-Lab/MetaScreener/tree/main/data
# Load DiffMetaScreener example data
load("diffGSE33113.Rdata")
load("diffGSE33113_G.Rdata")
load("diffGSE35896.Rdata")
load("diffGSE35896_G.Rdata")
load("diffGSE13067.Rdata")
load("diffGSE13067_G.Rdata")
load("example_diffgroup.Rdata")
load("example_geneset.Rdata")

# Load CorMetaScreener example data
load("corGSE33113.Rdata")
load("corGSE33113_V.Rdata")
load("corGSE35896.Rdata")
load("corGSE35896_V.Rdata")
load("corGSE13067.Rdata")
load("corGSE13067_V.Rdata")
```

## Usage Examples

### Mode 1: DiffMetaScreener - Based on Differential Analysis

This mode is designed for analyzing group comparisons (e.g., case vs. control, high vs. low expression) across multiple datasets.

```r
result_WNT_DiffMetaScreener <- DiffMetaScreener(
  # Vector of dataset names (must match loaded objects)
  expression_accession_vector = c("diffGSE33113", "diffGSE35896", "diffGSE13067"),
  
  # Comparison file specifying group labels and biological status
  comparisons = example_diffgroup,
  
  # Gene sets to be screened (as a GeneSetCollection object)
  geneSets_gmt = example_geneset,
  
  # Enrichment methods to apply
  enrichment_method = c(
    "gsva_t", "gsva_limma", "gsva_wilcoxon",
    "ssgsea_t", "ssgsea_limma", "ssgsea_wilcoxon",
    "zscore_t", "zscore_limma", "zscore_wilcoxon"
  ),
  
  # Minimum and maximum gene set sizes
  min.sz = 2,
  max.sz = 10000,
  
  # P-value combination methods
  p_combine_method = c("meanp", "meanz", "sumlog", "sumz", "sump")
)
```

### Mode 2: CorMetaScreener - Based on Correlation Analysis

This mode is designed for analyzing relationships between gene signatures and continuous variables (e.g., clinical scores, drug response) across multiple datasets.

```r
result_WNT_CorMetaScreener <- CorMetaScreener(
  # Vector of dataset names (must match loaded objects)
  expression_accession_vector = c("corGSE33113", "corGSE35896", "corGSE13067"),
  
  # Gene sets to be screened
  geneSets_gmt = example_geneset,
  
  # Correlation-based enrichment methods
  enrichment_method = c(
    "gsva_pearson", "gsva_kendall", "gsva_spearman",
    "ssgsea_pearson", "ssgsea_kendall", "ssgsea_spearman",
    "zscore_pearson", "zscore_kendall", "zscore_spearman"
  ),
  
  # Minimum and maximum gene set sizes
  min.sz = 2,
  max.sz = 10000,
  
  # P-value combination methods
  p_combine_method = c("meanp", "meanz", "sumlog", "sumz", "sump")
)
```

## Key Features

### 1. **Multiple Enrichment Methods**

MetaScreener supports a wide range of gene set enrichment methods:

- **GSVA** (Gene Set Variation Analysis)
- **ssGSEA** (single-sample GSEA)
- **z-score** methods
- **PLAGE** (Pathway Level Analysis of Gene Expression)
- **PCA**-based methods
- **AUCell** and **UCell**
- **singscore**
- **FGSEA** and **ORA** (Over-Representation Analysis)
- **decoupleR** methods (consensus, MLM, UDT, ULM, VIPER, WMEAN, WSUM)

### 2. **Statistical Approaches**

Each enrichment method can be combined with various statistical tests:

- **For DiffMetaScreener**: t-test, limma, Wilcoxon test, ANOVA, permutation test, Kruskal-Wallis
- **For CorMetaScreener**: Pearson, Kendall, Spearman correlation, linear regression, biweight, distance correlation, and other robust correlation methods

### 3. **P-value Integration**

Multiple methods for combining p-values across datasets and methods:

- Geometric mean
- Inverse chi-square (invchisq)
- Inverse t (invt)
- Logit method (logitp)
- Cauchy combination test (cct)
- Mean p-value (meanp)
- Mean z-score (meanz)
- Sum of logs (sumlog)
- Sum of z-scores (sumz)
- Sum of p-values (sump)
- Vote counting (votep)
- Wilkinson's method (wilkinsonp)

### 4. **Output Metrics**

The function returns two key metrics for each gene signature:

- **ADI (Activation Direction Index)**: Measures the strength of evidence for signature activation
- **IDI (Inhibition Direction Index)**: Measures the strength of evidence for signature inhibition

## Output Format

The result is a data frame with three columns:

```r
> head(result_WNT_DiffMetaScreener)
    signatures       ADI        IDI
1 WNT_signature1  2.345678  -1.234567
2 WNT_signature2  1.876543  -0.987654
3 ...            ...        ...
```

- **signatures**: Names of the gene sets
- **ADI**: Positive values indicate activation tendency
- **IDI**: Negative values indicate inhibition tendency

Higher absolute values indicate stronger and more consistent directional signals across datasets and methods.

## Citation

If you use MetaScreener in your research, please cite:

Zhao, D., Zhao, G., Yao, M. et al. MetaScreener: a robust dual-mode framework for directional prioritization of actionable signatures through multi-dataset and multi-approach integration. J Transl Med (2026). https://doi.org/10.1186/s12967-026-08019-y

## Author

Should you have any queries regarding installation or usage, please contact us by email. During spare moments, we can assist in calculating the directional scores for signatures.
Dingkang Zhao (dingkang.25@intl.zju.edu.cn)

## License

This package is distributed under the terms specified in the package description.
