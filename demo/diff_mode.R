# DiffMetaScreener(): directional screening against a two-group phenotype
#
# MetaScreener scores each signature in every cohort with a panel of enrichment
# methods, tests the scores between the activation and inhibition groups,
# converts the two-tailed p-values into directional one-sided evidence and
# combines that evidence across methods and datasets.
#
# This demo uses two of the 23 enrichment methods so that it finishes in about
# twenty seconds. Drop the `enrichment_method` argument to run the full set.

library(MetaScreener)

data("example_geneset")     # 10 Wnt signatures
data("example_diffgroup")   # which groups to compare in which dataset
data(list = c("diffGSE33113", "diffGSE33113_G",
              "diffGSE35896", "diffGSE35896_G",
              "diffGSE13067", "diffGSE13067_G"))

cat("\nThree cohorts:", paste(dim(diffGSE33113)[1], "genes x",
                              dim(diffGSE33113)[2], "samples\n"))
cat("Groups:", paste(names(table(diffGSE33113_G$group)),
                     table(diffGSE33113_G$group), sep = "=", collapse = ", "), "\n\n")

res <- DiffMetaScreener(
  expression_accession_vector = c("diffGSE33113", "diffGSE35896", "diffGSE13067"),
  comparisons                 = example_diffgroup,
  geneSets_gmt                = example_geneset,
  enrichment_method           = c("zscore_t", "zscore_wilcoxon"),
  min.sz = 2,
  max.sz = 10000
)

res$DI <- res$ADI + res$IDI
print(res[order(-res$DI), ])

cat("\nADI = mean -log10 of the combined activation p-values (always positive)",
    "\nIDI = mean  log10 of the combined inhibition p-values (always negative)",
    "\nDI  = ADI + IDI, the net directionality used for ranking\n")
