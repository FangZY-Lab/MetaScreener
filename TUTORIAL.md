# MetaScreener — a step-by-step tutorial

This tutorial takes you from an empty R session to a ranked table of directional
signatures. It assumes you have never used the package before. Every step shows
what to type and what you should see.

Total time: about 10 minutes, most of it spent installing dependencies once.

---

## Before you start

You need:

- R ≥ 4.0 (4.3 or newer recommended);
- RStudio or any editor that can run R;
- an internet connection for the one-off installation.

The package pulls in a fairly large stack of Bioconductor packages (GSVA,
clusterProfiler, AUCell, UCell, singscore, decoupleR, fgsea, limma, GSEABase),
plus `IOBR` from GitHub and a handful of CRAN packages. The first installation
takes 15–40 minutes depending on your machine. Later steps are fast.

---

## Step 0 — Install once

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
install.packages("remotes")

# IOBR lives on GitHub, not on CRAN or Bioconductor
remotes::install_github("IOBR/IOBR")

BiocManager::install(c(
  "GSVA", "GSEABase", "AUCell", "singscore", "UCell",
  "decoupleR", "clusterProfiler", "limma", "fgsea"
), update = FALSE)

install.packages(c("metap", "doBy", "dplyr", "plyr", "reshape2"))

remotes::install_github("FangZY-Lab/MetaScreener")
```

Check that it worked:

```r
library(MetaScreener)
packageVersion("MetaScreener")     # should print 1.1.0 or later
```

> **If installation fails** with `there is no package called 'X'`, see
> [Troubleshooting](#troubleshooting) at the bottom.

---

## Step 1 — Load the example data

The package ships with three colorectal cancer cohorts, two complete example
sets and a `GeneSetCollection` of ten Wnt signatures. Load them with `data()`:

```r
data("example_geneset")     # 10 signatures
data("example_diffgroup")   # which groups to compare in which dataset

data(list = c(                              # the three Diff-mode cohorts
  "diffGSE33113", "diffGSE33113_G",
  "diffGSE35896", "diffGSE35896_G",
  "diffGSE13067", "diffGSE13067_G"
))
```

`data()` copies the objects into your global environment, which is where the
screener looks for them. If you ever get `object 'diffGSE33113' not found`,
it is because this step was skipped.

**Check the objects arrived:**

```r
dim(diffGSE33113)        # 21502  80
head(diffGSE33113_G, 3)  # Tag, group
```

You should see a matrix with 21,502 genes in the rows and 80 samples in the
columns, and a two-column table whose `Tag` column matches the matrix columns.

---

## Step 2 — Look at what you are screening

```r
example_diffgroup
#>      Accession Group1 Group2 Group1_Status Group2_Status
#> 1 diffGSE33113   CMS2  Other    activation    inhibition
#> 2 diffGSE35896   CMS2  Other    activation    inhibition
#> 3 diffGSE13067   CMS2  Other    activation    inhibition

length(example_geneset)          # 10 signatures
head(names(example_geneset))     # BIOCARTA_WNT ...
```

`comparisons` tells the screener, for each dataset, which two groups to compare
and which one represents **activation**. In the example, CMS2 is the activation
group and everything else is the inhibition group.

---

## Step 3 — Run your first screen

Start small: three datasets, two cheap scoring methods, default p-value
combination. This takes about 20 seconds.

```r
res <- DiffMetaScreener(
  expression_accession_vector = c("diffGSE33113", "diffGSE35896", "diffGSE13067"),
  comparisons                 = example_diffgroup,
  geneSets_gmt                = example_geneset,
  enrichment_method           = c("zscore_t", "zscore_wilcoxon"),
  min.sz = 2,
  max.sz = 10000
)
```

The console prints a lot of progress (`gsva_t will start.`, `Not using viper.`,
…). That is normal: the function reports every method it tries.

---

## Step 4 — Read the output

```r
res
#>      signatures         ADI           IDI
#> 1  BIOCARTA_WNT 12.01998509 -2.179156e-02
#> 2   FEVR_WNT_DN  0.01452660 -8.911991e+00
#> 3   FEVR_WNT_UP  5.30784588 -3.516706e-02
#> 4     FLIER_WNT 18.74944260 -6.136770e-06
#> 5      GOBP_WNT  0.57010137 -1.005907e+00
#> ...
```

Three columns:

| Column | Meaning |
| --- | --- |
| `signatures` | the name of the gene set |
| `ADI` | **A**ctivation **D**irection **I**ndex — mean of −log₁₀ of the combined activation p-values, so always positive |
| `IDI` | **I**nhibition **D**irection **I**ndex — mean of log₁₀ of the combined inhibition p-values, so always negative |

The paper also defines the net directionality index, `DI = ADI + IDI`, which is
what ranking should use. The functions return `ADI` and `IDI`; add the column
yourself:

```r
res$DI <- res$ADI + res$IDI
```

How to read a row:

- **large ADI, IDI near zero** → the signature is consistently activated.
  `FLIER_WNT` (ADI 18.7) is the clearest example.
- **ADI near zero, large negative IDI** → consistently inhibited.
  `FEVR_WNT_DN` (IDI −8.9) is the clearest example.
- **ADI and IDI both substantial, opposite signs** → context-dependent; the sign
  of `DI` tells you which direction wins overall.
- **both near zero** → no consistent direction; drop it or look again with more
  datasets.

To rank:

```r
res[order(-res$DI), ]                     # net activation first
res[order(-res$ADI), ]                    # strongest activation evidence
res[order(res$IDI), ]                     # strongest inhibition evidence
```

---

## Step 5 — Run the correlation mode

`CorMetaScreener()` answers a different question. Instead of two labelled groups
it takes a **continuous state variable** per sample — a pathway activity score, a
clinical measurement, a drug sensitivity value — and asks whether a signature's
activity tracks it. The bundled `cor` objects carry such a variable in their
`_V` tables.

```r
data(list = c("corGSE33113", "corGSE33113_V",
              "corGSE35896", "corGSE35896_V",
              "corGSE13067", "corGSE13067_V"))

res_cor <- CorMetaScreener(
  expression_accession_vector = c("corGSE33113", "corGSE35896", "corGSE13067"),
  geneSets_gmt                = example_geneset,
  enrichment_method           = c("zscore_pearson", "zscore_spearman"),
  min.sz = 2,
  max.sz = 10000
)

res_cor$DI <- res_cor$ADI + res_cor$IDI
res_cor[order(-res_cor$DI), ]
```

```text
   signatures               ADI            IDI             DI
1  HALLMARK_WNT    16.710592822  -2.774547e-05   16.710565077
2  WP_WNT          13.928138294  -1.593803e-03   13.926544491
3  FLIER_WNT       12.749898469  -6.007688e-08   12.749898409
4  BIOCARTA_WNT     5.796530160  -4.734137e-03    5.791796022
5  REACTOME_WNT     5.579781384  -1.292526e-01    5.450528833
6  KEGG_WNT         4.844878924  -3.839505e-02    4.806483874
7  FEVR_WNT_UP      1.638466082  -7.592722e-02    1.562538860
8  GOBP_WNT         0.736208787  -1.711992e-01    0.565009617
9  PID_WNT          1.258082003  -3.555518e+00   -2.297436086
10 FEVR_WNT_DN      0.003118185  -1.146383e+01  -11.460715153
```

Two things are worth noticing:

- the output has exactly the same three columns as the differential run, so
  everything in Step 4 applies unchanged;
- the ranking is not the same. `HALLMARK_WNT` and `WP_WNT` lead here, `FLIER_WNT`
  led the differential screen, and `FEVR_WNT_DN` is inhibited under both. That is
  expected: the two modes test different hypotheses on the same signatures, and
  reporting both is stronger than reporting either alone.

## Step 6 — Screen your own data

You need one expression matrix and one annotation table per dataset. Everything
is matched by sample ID, so the annotation table must have a `Tag` column whose
values equal `colnames(<matrix>)`.

```r
# a minimal, correct dataset
my_expr <- as.matrix(read.csv("my_cohort_expr.csv", row.names = 1))   # genes x samples
my_ann  <- read.csv("my_cohort_annotation.csv")                        # columns Tag, group
#   Tag   group
#   S1    CMS2
#   S2    Other

stopifnot(setequal(my_ann$Tag, colnames(my_expr)))   # alignment check

# the screener looks objects up by name, so name them <id> and <id>_G
myCohort      <- my_expr
myCohort_G    <- my_ann

# and describe the comparison
my_comparisons <- data.frame(
  Accession     = "myCohort",
  Group1        = "CMS2",
  Group2        = "Other",
  Group1_Status = "activation",
  Group2_Status = "inhibition"
)

res_mine <- DiffMetaScreener(
  expression_accession_vector = "myCohort",
  comparisons                 = my_comparisons,
  geneSets_gmt                = example_geneset,   # or your own GeneSetCollection
  enrichment_method           = c("zscore_t", "zscore_wilcoxon"),
  min.sz = 2, max.sz = 10000
)
```

For the correlation mode the annotation table has columns `Tag` and `value`
instead, and the object is named `<id>_V`:

```r
myValue   <- data.frame(Tag = my_ann$Tag, value = my_ann$score)
corMyCohort   <- my_expr
corMyCohort_V <- myValue

res_cor <- CorMetaScreener(
  expression_accession_vector = "corMyCohort",
  geneSets_gmt                = example_geneset,
  enrichment_method           = c("zscore_pearson", "zscore_spearman"),
  min.sz = 2, max.sz = 10000
)
```

Load a GMT file into a `GeneSetCollection` with:

```r
my_sets <- GSEABase::getGmt("my_signatures.gmt")
```

---

## Step 7 — Choose your methods

Two arguments control the work:

**`enrichment_method`** — how each signature's activity is scored. The name has
two halves separated by `_`: the scoring family and the statistical test.

```text
zscore_t              z-score scoring, t-test
gsva_wilcoxon         GSVA scoring, Wilcoxon test
singscore_limma       singscore, limma
aucell_permutation    AUCell, permutation test
```

Scoring families: `gsva`, `ssgsea`, `zscore`, `plage`, `pca`, `aucell`, `ucell`,
`singscore`, `median`, `fgsea`, `ora`, `consensus`, `mdt`, `mlm`, `udt`, `ulm`,
`viper`, `wmean`, `norm_wmean`, `corr_wmean`, `wsum`, `norm_wsum`, `corr_wsum`.

Tests for `DiffMetaScreener`: `t`, `limma`, `anova`, `wilcoxon`, `permutation`,
`kruskal`. For `CorMetaScreener`: `pearson`, `spearman`, `kendall`, `lm`,
`biweight`, `distance`, `percentage`, `blomqvist`, `hoeffding`, `gamma`.

**`p_combine_method`** — how p-values are merged across the enrichment methods
**within each dataset** (twelve meta-analysis algorithms from `metap` plus CCT
and the geometric mean). The cross-dataset step is fixed: it always uses the
geometric mean. Leave this at the default (all twelve).

> **Do not** narrow it to `c("meanp", "meanz")` with few methods and few
> datasets. Both need at least four valid p-values at *each* combining step;
> below that they return `NA`, the pipeline substitutes 1, and every ADI/IDI
> becomes exactly 0. The function now warns you when this is about to happen.

---

## Step 8 — Scale up

The full default set of methods is large. A practical progression:

| Stage | Methods | Datasets | Rough runtime |
| --- | --- | --- | --- |
| first look | 2 | 2–3 | under a minute |
| working screen | 6–10 | 3–5 | several minutes |
| publication run | all defaults | all cohorts | hours |

Runtime grows linearly with datasets × methods, so test the pipeline on a small
subset before launching the full job.

---

## Troubleshooting

| Message | Cause | Fix |
| --- | --- | --- |
| `object 'diffGSE33113' not found` | `data()` not called | run the `data()` calls from Step 1 |
| `data set 'diffGSE33113' not found` | old package version, or data files missing | reinstall with `remotes::install_github("FangZY-Lab/MetaScreener")` |
| `Missing objects for: X. Each dataset needs ...` | object name typo | the matrix must be `<id>`, the annotation `<id>_G` or `<id>_V` |
| `` `X`$Tag and colnames(`X`) do not contain the same samples `` | annotation and matrix disagree | fix the IDs so both sets match exactly |
| `there is no package called 'X'` | dependency missing | `BiocManager::install("X")` or `install.packages("X")` |
| every ADI/IDI is exactly 0 | `meanp`/`meanz` with fewer than 4 p-values | keep the default `p_combine_method` |
| `Must have at least four valid p values` warning | same cause | same fix |
| results are all `NA` | a dataset has too few samples per group | check `table(dataset_G$group)` |

Getting help: open an issue at
<https://github.com/FangZY-Lab/MetaScreener/issues> with the output of
`sessionInfo()` and the call you ran.
