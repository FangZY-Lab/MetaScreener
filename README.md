# MetaScreener

## Overview

**MetaScreener: a robust dual-mode framework for directional prioritization of actionable signatures through multi-dataset and multi-approach integration**

MetaScreener represents a paradigm shift in gene set enrichment analysis by introducing a sophisticated multi-layer meta-analysis framework that transcends traditional single-dataset, single-method approaches. The package's groundbreaking architecture enables:

**Core Innovation: DiffMetaScreener**
A revolutionary differential enrichment analysis engine that performs directional signature screening through:
- **Multi-method agnostic integration**: Simultaneously employs 60+ enrichment methodologies spanning gene set scoring, network inference, and pathway enrichment approaches
- **Directional p-value transformation**: Converts conventional two-tailed p-values into activation-specific and inhibition-specific evidence metrics
- **Hierarchical meta-analysis**: Executes two-tiered evidence integration across both methodological and dataset dimensions
- **Dual directional indices**: Generates Activation Direction Index (ADI) and Inhibition Direction Index (IDI) for unambiguous biological interpretation

**Secondary Mode: CorMetaScreener**
A correlation-based enrichment module that extends the multi-layer meta-analysis principle to continuous phenotypic associations.

**Transformative Capabilities:**
- Robust cross-study signature validation through methodological consensus
- Elimination of single-method bias via evidence aggregation
- Directional pathway activity quantification beyond conventional differential expression
- Scalable framework accommodating diverse experimental designs and technological platforms

MetaScreener moves beyond conventional enrichment analysis by providing a systematic, evidence-based framework for discovering and validating directional gene set signatures across complex multi-omics landscapes.

## Installation

### Prerequisites

Make sure you have R (version 4.0.0 or higher) installed. The package requires several Bioconductor and CRAN dependencies.

### Installation Steps

```r
# Clear workspace and set working directory
rm(list=ls())
gc()
setwd("your/working/directory/path") # Set your working directory here

# Install dependency packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# List of required packages
depens = c("singscore","UCell","decoupleR","AUCell","metap","IOBR",
           "fgsea","plyr","reshape2","GSEABase","GSVA","doBy",
           "clusterProfiler","limma","dplyr","tidyverse","devtools","stringr",
           "GSA","coin","viper","caret","wdm","correlation")

# Install missing dependencies
for(i in 1:length(depens)){
  depen = depens[i]
  if (!requireNamespace(depen, quietly = TRUE))
    BiocManager::install(depen, update = FALSE)
}

# Install MetaScreener from GitHub
devtools::install_github("FangZY-Lab/MetaScreener", force = TRUE)

# Load the package
library(MetaScreener)
# Load example data for DiffMetaScreener mode
# The data is available for download in the 'data' directory of the FangZY-Lab/MetaScreener repository on GitHub
load("example_diffdata.Rdata")
load("example_diffdata_G.Rdata")
load("example_diffgroup.Rdata")
load("example_geneset.Rdata")

# Load example data for CorMetaScreener mode
# The data is available for download in the 'data' directory of the FangZY-Lab/MetaScreener repository on GitHub
load("example_cordata.Rdata")
load("example_cordata_V.Rdata")
# Assign loaded data to variables
for(i in 1:length(names(example_diffdata))){
  exp = example_diffdata[[i]]
  assign(names(example_diffdata)[i], exp)
}
rm(exp)

for(i in 1:length(names(example_diffdata_G))){
  exp_G = example_diffdata_G[[i]]
  assign(names(example_diffdata_G)[i], exp_G)
}
rm(exp_G)

# Run DiffMetaScreener analysis
result_WNT_DiffMetaScreener = DiffMetaScreener(
  expression_accession_vector = names(example_diffdata),
  comparisons = example_diffgroup,
  geneSets_gmt = example_geneset,
  enrichment_method = c(
    "gsva_t","gsva_limma","gsva_anova","gsva_wilcoxon","gsva_permutation","gsva_kruskal",
    "ssgsea_t","ssgsea_limma","ssgsea_anova","ssgsea_wilcoxon","ssgsea_permutation","ssgsea_kruskal",
    "zscore_t","zscore_limma","zscore_anova","zscore_wilcoxon","zscore_permutation","zscore_kruskal",
    "plage_t","plage_limma","plage_anova","plage_wilcoxon","plage_permutation","plage_kruskal",
    "pca_t","pca_limma","pca_anova","pca_wilcoxon","pca_permutation","pca_kruskal",
    "aucell_t","aucell_limma","aucell_anova","aucell_wilcoxon","aucell_permutation","aucell_kruskal",
    "ucell_t","ucell_limma","ucell_anova","ucell_wilcoxon","ucell_permutation","ucell_kruskal",
    "singscore_t","singscore_limma","singscore_anova","singscore_wilcoxon","singscore_permutation","singscore_kruskal",
    "median_t","median_limma","median_anova","median_wilcoxon","median_permutation","median_kruskal",
    "t_fgsea","limma_fgsea","anova_fgsea","wilcoxon_fgsea","permutation_fgsea","kruskal_fgsea",
    "t_ora","limma_ora","anova_ora","wilcoxon_ora","permutation_ora","kruskal_ora",
    "consensus_t","consensus_limma","consensus_anova","consensus_wilcoxon","consensus_permutation","consensus_kruskal",
    "mdt_t","mdt_limma","mdt_anova","mdt_wilcoxon","mdt_permutation","mdt_kruskal",
    "mlm_t","mlm_limma","mlm_anova","mlm_wilcoxon","mlm_permutation","mlm_kruskal",
    "udt_t","udt_limma","udt_anova","udt_wilcoxon","udt_permutation","udt_kruskal",
    "ulm_t","ulm_limma","ulm_anova","ulm_wilcoxon","ulm_permutation","ulm_kruskal",
    "viper_t","viper_limma","viper_anova","viper_wilcoxon","viper_permutation","viper_kruskal",
    "wmean_t","wmean_limma","wmean_anova","wmean_wilcoxon","wmean_permutation","wmean_kruskal",
    "norm_wmean_t","norm_wmean_limma","norm_wmean_anova","norm_wmean_wilcoxon","norm_wmean_permutation","norm_wmean_kruskal",
    "corr_wmean_t","corr_wmean_limma","corr_wmean_anova","corr_wmean_wilcoxon","corr_wmean_permutation","corr_wmean_kruskal",
    "wsum_t","wsum_limma","wsum_anova","wsum_wilcoxon","wsum_permutation","wsum_kruskal",
    "norm_wsum_t","norm_wsum_limma","norm_wsum_anova","norm_wsum_wilcoxon","norm_wsum_permutation","norm_wsum_kruskal",
    "corr_wsum_t","corr_wsum_limma","corr_wsum_anova","corr_wsum_wilcoxon","corr_wsum_permutation","corr_wsum_kruskal"
  ),
  min.sz = 2,
  max.sz = 10000,
  p_combine_method = c(
    "geometric_mean","invt","invchisq","logitp","cct",
    "meanp","meanz","sumlog","sumz","sump","votep","wilkinsonp"
  )
)
# Assign loaded data to variables
for(i in 1:length(names(example_cordata))){
  exp = example_cordata[[i]]
  assign(names(example_cordata)[i], exp)
}
rm(exp)

for(i in 1:length(names(example_cordata_V))){
  exp_V = example_cordata_V[[i]]
  assign(names(example_cordata_V)[i], exp_V)
}
rm(exp_V)

# Run CorMetaScreener analysis
result_WNT_CorMetaScreener = CorMetaScreener(
  expression_accession_vector = names(example_cordata),
  geneSets_gmt = example_geneset,
  enrichment_method = c(
    "gsva_pearson","gsva_kendall","gsva_spearman","gsva_lm","gsva_biweight","gsva_distance","gsva_percentage","gsva_blomqvist","gsva_hoeffding","gsva_gamma",
    "ssgsea_pearson","ssgsea_kendall","ssgsea_spearman","ssgsea_lm","ssgsea_biweight","ssgsea_distance","ssgsea_percentage","ssgsea_blomqvist","ssgsea_hoeffding","ssgsea_gamma",
    "zscore_pearson","zscore_kendall","zscore_spearman","zscore_lm","zscore_biweight","zscore_distance","zscore_percentage","zscore_blomqvist","zscore_hoeffding","zscore_gamma",
    "plage_pearson","plage_kendall","plage_spearman","plage_lm","plage_biweight","plage_distance","plage_percentage","plage_blomqvist","plage_hoeffding","plage_gamma",
    "pca_pearson","pca_kendall","pca_spearman","pca_lm","pca_biweight","pca_distance","pca_percentage","pca_blomqvist","pca_hoeffding","pca_gamma",
    "aucell_pearson","aucell_kendall","aucell_spearman","aucell_lm","aucell_biweight","aucell_distance","aucell_percentage","aucell_blomqvist","aucell_hoeffding","aucell_gamma",
    "ucell_pearson","ucell_kendall","ucell_spearman","ucell_lm","ucell_biweight","ucell_distance","ucell_percentage","ucell_blomqvist","ucell_hoeffding","ucell_gamma",
    "singscore_pearson","singscore_kendall","singscore_spearman","singscore_lm","singscore_biweight","singscore_distance","singscore_percentage","singscore_blomqvist","singscore_hoeffding","singscore_gamma",
    "median_pearson","median_kendall","median_spearman","median_lm","median_biweight","median_distance","median_percentage","median_blomqvist","median_hoeffding","median_gamma",
    "pearson_fgsea","kendall_fgsea","spearman_fgsea","lm_fgsea","biweight_fgsea","distance_fgsea","percentage_fgsea","blomqvist_fgsea","hoeffding_fgsea","gamma_fgsea",
    "pearson_ora","kendall_ora","spearman_ora","lm_ora","biweight_ora","distance_ora","percentage_ora","blomqvist_ora","hoeffding_ora","gamma_ora",
    "consensus_pearson","consensus_kendall","consensus_spearman","consensus_lm","consensus_biweight","consensus_distance","consensus_percentage","consensus_blomqvist","consensus_hoeffding","consensus_gamma",
    "mdt_pearson","mdt_kendall","mdt_spearman","mdt_lm","mdt_biweight","mdt_distance","mdt_percentage","mdt_blomqvist","mdt_hoeffding","mdt_gamma",
    "mlm_pearson","mlm_kendall","mlm_spearman","mlm_lm","mlm_biweight","mlm_distance","mlm_percentage","mlm_blomqvist","mlm_hoeffding","mlm_gamma",
    "udt_pearson","udt_kendall","udt_spearman","udt_lm","udt_biweight","udt_distance","udt_percentage","udt_blomqvist","udt_hoeffding","udt_gamma",
    "ulm_pearson","ulm_kendall","ulm_spearman","ulm_lm","ulm_biweight","ulm_distance","ulm_percentage","ulm_blomqvist","ulm_hoeffding","ulm_gamma",
    "viper_pearson","viper_kendall","viper_spearman","viper_lm","viper_biweight","viper_distance","viper_percentage","viper_blomqvist","viper_hoeffding","viper_gamma",
    "wmean_pearson","wmean_kendall","wmean_spearman","wmean_lm","wmean_biweight","wmean_distance","wmean_percentage","wmean_blomqvist","wmean_hoeffding","wmean_gamma",
    "norm_wmean_pearson","norm_wmean_kendall","norm_wmean_spearman","norm_wmean_lm","norm_wmean_biweight","norm_wmean_distance","norm_wmean_percentage","norm_wmean_blomqvist","norm_wmean_hoeffding","norm_wmean_gamma",
    "corr_wmean_pearson","corr_wmean_kendall","corr_wmean_spearman","corr_wmean_lm","corr_wmean_biweight","corr_wmean_distance","corr_wmean_percentage","corr_wmean_blomqvist","corr_wmean_hoeffding","corr_wmean_gamma",
    "wsum_pearson","wsum_kendall","wsum_spearman","wsum_lm","wsum_biweight","wsum_distance","wsum_percentage","wsum_blomqvist","wsum_hoeffding","wsum_gamma",
    "norm_wsum_pearson","norm_wsum_kendall","norm_wsum_spearman","norm_wsum_lm","norm_wsum_biweight","norm_wsum_distance","norm_wsum_percentage","norm_wsum_blomqvist","norm_wsum_hoeffding","norm_wsum_gamma",
    "corr_wsum_pearson","corr_wsum_kendall","corr_wsum_spearman","corr_wsum_lm","corr_wsum_biweight","corr_wsum_distance","corr_wsum_percentage","corr_wsum_blomqvist","corr_wsum_hoeffding","corr_wsum_gamma"
  ),
  min.sz = 2,
  max.sz = 10000,
  p_combine_method = c(
    "geometric_mean","invt","invchisq","logitp","cct",
    "meanp","meanz","sumlog","sumz","sump","votep","wilkinsonp"
  )
)
