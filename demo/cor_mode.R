# CorMetaScreener(): directional screening against a continuous state variable
#
# The correlative mode asks whether a signature's activity tracks a continuous
# variable attached to each sample instead of separating two labelled groups.
#
# This demo uses two of the 23 enrichment methods so that it finishes in about
# twenty seconds. Drop the `enrichment_method` argument to run the full set.

library(MetaScreener)

data("example_geneset")     # 10 Wnt signatures
data(list = c("corGSE33113", "corGSE33113_V",
              "corGSE35896", "corGSE35896_V",
              "corGSE13067", "corGSE13067_V"))

cat("\nThree cohorts:", paste(dim(corGSE33113)[1], "genes x",
                              dim(corGSE33113)[2], "samples\n"))
cat("State variable range:",
    paste(round(range(corGSE33113_V$value), 3), collapse = " to "), "\n\n")

res_cor <- CorMetaScreener(
  expression_accession_vector = c("corGSE33113", "corGSE35896", "corGSE13067"),
  geneSets_gmt                = example_geneset,
  enrichment_method           = c("zscore_pearson", "zscore_spearman"),
  min.sz = 2,
  max.sz = 10000
)

res_cor$DI <- res_cor$ADI + res_cor$IDI
print(res_cor[order(-res_cor$DI), ])

cat("\nNote that the ranking differs from demo('diff_mode'): the two modes test",
    "\ndifferent hypotheses on the same signatures.\n")
