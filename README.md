# MetaScreener

**Directional, multi-dataset screening of gene signatures.**

MetaScreener takes a set of gene signatures and asks a single question of each
one: across many independent datasets and many enrichment methods, is the
signature consistently *activated* or consistently *inhibited*? It answers with
two scores per signature — an **Activation Direction Index (ADI)** and an
**Inhibition Direction Index (IDI)** — rather than a single p-value that hides
the direction of the effect.

Two modes share the same engine:

| Mode | Input | Typical use |
| --- | --- | --- |
| `DiffMetaScreener()` | expression + a two-group annotation per dataset | case vs. control, subtype A vs. subtype B, high vs. low |
| `CorMetaScreener()` | expression + a continuous variable per dataset | clinical score, drug sensitivity, pathway activity |

## Installation

MetaScreener depends on several Bioconductor packages, and on `IOBR`, which is
distributed through GitHub. Install those first:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
install.packages("remotes")

# IOBR is not on CRAN or Bioconductor
remotes::install_github("IOBR/IOBR")

BiocManager::install(c(
  "GSVA", "GSEABase", "AUCell", "singscore", "UCell",
  "decoupleR", "clusterProfiler", "limma", "fgsea"
), update = FALSE)

install.packages(c("metap", "doBy", "dplyr", "plyr", "reshape2"))
```

Then install the package itself:

```r
remotes::install_github("FangZY-Lab/MetaScreener")
library(MetaScreener)
```

## Quick start

The package ships with two complete example sets (three datasets each) and a
ten-signature `GeneSetCollection`, so the workflow below runs as written.

```r
library(MetaScreener)

# signatures and the comparison table
data("example_geneset")     # GeneSetCollection, 10 signatures
data("example_diffgroup")   # which groups to compare in each dataset

# Diff mode: three colorectal cohorts, CMS2 versus the rest
data(list = c(
  "diffGSE33113", "diffGSE33113_G",
  "diffGSE35896", "diffGSE35896_G",
  "diffGSE13067", "diffGSE13067_G"
))
```

`data()` loads the objects into the global environment, which is where the
screener looks for them.

### Mode 1 — differential analysis

```r
result_diff <- DiffMetaScreener(
  expression_accession_vector = c("diffGSE33113", "diffGSE35896", "diffGSE13067"),
  comparisons                 = example_diffgroup,
  geneSets_gmt                = example_geneset,
  enrichment_method = c("gsva_t", "gsva_wilcoxon",
                        "ssgsea_t", "ssgsea_wilcoxon",
                        "zscore_t", "zscore_wilcoxon"),
  min.sz = 2,
  max.sz = 10000,
  p_combine_method = c("meanp", "meanz", "sumlog", "sumz", "sump")
)

head(result_diff)
```

### Mode 2 — correlation analysis

```r
data(list = c(
  "corGSE33113", "corGSE33113_V",
  "corGSE35896", "corGSE35896_V",
  "corGSE13067", "corGSE13067_V"
))

result_cor <- CorMetaScreener(
  expression_accession_vector = c("corGSE33113", "corGSE35896", "corGSE13067"),
  geneSets_gmt = example_geneset,
  enrichment_method = c("gsva_pearson", "gsva_spearman",
                        "ssgsea_pearson", "ssgsea_spearman",
                        "zscore_pearson", "zscore_spearman"),
  min.sz = 2,
  max.sz = 10000,
  p_combine_method = c("meanp", "meanz", "sumlog", "sumz", "sump")
)

head(result_cor)
```

Both functions accept around a hundred `enrichment_method` values; the short
lists above keep a first run quick. Run the full default set only when you are
ready for a long job — every method is fitted in every dataset.

## Input contract

The screener reads objects from the **global environment**, by name. For each
entry in `expression_accession_vector` the following must exist:

| Object | Content |
| --- | --- |
| `<id>` | expression matrix, genes in rows, samples in columns, log2 scale |
| `<id>_G` (Diff) | data frame with columns `Tag` (sample ID) and `group` |
| `<id>_V` (Cor) | data frame with columns `Tag` (sample ID) and `value` |

`Tag` must match `colnames(<id>)` exactly — same values, and the bundled data
keeps the same order. Row names of the expression matrix must be gene symbols
without duplicates, and the gene sets must use the same identifier type.

`comparisons` is a data frame with one row per dataset:

| Accession | Group1 | Group2 | Group1_Status | Group2_Status |
| --- | --- | --- | --- | --- |
| diffGSE33113 | CMS2 | Other | activation | inhibition |

`GeneSetCollection` objects can be read with `GSEABase::getGmt()` from a GMT
file, which is the usual starting point for your own signatures.

## Output

Both functions return a data frame with one row per signature:

```text
  signatures        ADI        IDI
1 WNT_signature1   2.345678  -1.234567
2 WNT_signature2   1.876543  -0.987654
```

- `ADI` — evidence that the signature is activated (higher is stronger);
- `IDI` — evidence that it is inhibited;
- a large positive ADI with a small |IDI| means the signature behaves as an
  activation signature across datasets and methods, and vice versa.

## Bundled example data

| Object | Dimensions | Description |
| --- | --- | --- |
| `diffGSE33113`, `diffGSE35896`, `diffGSE13067` | 21502 genes x 80 / 54 / 58 samples | log2 expression, colorectal cancer |
| `<id>_G` | samples x 2 (`Tag`, `group`) | CMS2 or Other |
| `corGSE33113`, `corGSE35896`, `corGSE13067` | 21502 genes x 90 / 62 / 74 samples | log2 expression, colorectal cancer |
| `<id>_V` | samples x 2 (`Tag`, `value`) | continuous score per sample |
| `example_diffgroup` | 3 x 5 | comparison table for the three diff datasets |
| `example_geneset` | `GeneSetCollection` of 10 | signatures screened in the examples |

The three cohorts are public GEO series (GSE13067, GSE33113, GSE35896).

## Citation

Zhao D, Zhao G, Yao M, Wu J, Fang Z. MetaScreener: a robust dual-mode framework
for directional prioritization of actionable signatures through multi-dataset
and multi-approach integration. *J Transl Med*. 2026;24(1):583.
doi:10.1186/s12967-026-08019-y. PMID: 41845481.

## Contact

Dingkang Zhao — <dingkang.22@intl.zju.edu.cn>

## License

MIT. See `LICENSE`.
