# MetaScreener

**A dual-mode framework for directional prioritization of actionable signatures
through multi-dataset and multi-approach integration.**

[![License: MIT](https://img.shields.io/badge/license-MIT-1a7f37.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-%E2%89%A5%204.0-0969da.svg)](https://cran.r-project.org)

Reliably quantifying the directionality of gene signatures from noisy omics data
remains challenging. MetaScreener is the software implementation of the
framework described in [Zhao *et al.*, *J Transl Med* 2026](https://doi.org/10.1186/s12967-026-08019-y):
it integrates multiple datasets and a large number of analysis pipelines
(>4,000) to quantify **activating and inhibiting signals** in gene signatures,
and returns an activation directionality index (ADI) and an inhibition
directionality index (IDI) for every candidate signature.

![MetaScreener workflow](man/figures/MetaScreener-workflow.svg)

## The framework

MetaScreener was designed to satisfy four requirements: the compilation, quality
control, directional scoring and context-specific evaluation of candidate
signatures.

**Step 1 — Datasets and signatures.** One or more numerical matrices
(transcriptomics, epigenomics, proteomics, metabolomics — any calibrated omics
data whose features match the signatures). Each sample is annotated either with
binary *activation / inhibition* labels or with continuous state variables such
as gene expression, pathway activity, clinical features or user-defined values.
Candidate signatures may come from MSigDB, Enrichr, HMDB, PathCards, the
literature, custom definitions or de novo gene sets.

**Step 2 — Multi-pipeline aggregation.** Two modes share the same engine:

- `DiffMetaScreener()` — differential enrichment analysis on binary labels;
- `CorMetaScreener()` — correlative enrichment analysis on continuous variables.

Both combine **23 signature enrichment methods** — 21 unsupervised (GSVA,
ssGSEA, z-score, PLAGE, PCA, AUCell, UCell, singscore, median, consensus and the
eleven `decoupleR` approaches: MDT, MLM, UDT, ULM, VIPER, WMEAN, norm_WMEAN,
corr_WMEAN, WSUM, norm_WSUM, corr_WSUM) and 2 supervised (FGSEA, ranked by
−log₁₀*p* × log₂FC, and ORA on the top 1000 up- and down-regulated genes).

The enrichment scores are then tested per dataset:

| Mode | Tests | Pipelines |
| --- | --- | --- |
| Diff | limma, Student's *t*, ANOVA, Wilcoxon rank-sum, Kruskal–Wallis, permutation | 23 × 6 = **138** differential enrichment pipelines |
| Cor | Pearson, Kendall, Spearman, linear model, Biweight, Distance, Percentage, Blomqvist, Hoeffding, Gamma | 23 × 10 = **230** correlative enrichment pipelines |

That is **368 distinct functional profiling strategies** per analysis. The
two-tailed p-values are converted to directional one-sided p-values with
`two2one`, and combined **within each dataset** by 12 meta-analysis algorithms:
geometric mean, invt, invchisq, logitp, CCT, meanp, meanz, sumlog, sumz, sump,
votep and wilkinsonp. **Across datasets only the geometric mean is used.**

**Step 3 — Directionality indices.** The combined activation and inhibition
p-values are transformed by −log₁₀ and log₁₀ respectively, giving consistently
positive ADI and consistently negative IDI values for each of the twelve
algorithms. Averaging those values gives the final indices, and their sum is the
net directionality index used for ranking:

```text
ADI = mean over 12 algorithms of  -log10(p_activation)      always positive
IDI = mean over 12 algorithms of   log10(p_inhibition)      always negative
DI  = ADI + IDI                                             net directionality
```

The two functions return `signatures`, `ADI` and `IDI`; `DI` is their sum.

## What the paper reports

- **Reliability.** Ten candidate Wnt pathway signatures were screened in nine
  colorectal cancer training datasets. CMS2 samples (Wnt hyperactive) were the
  activation group and CMS1/3/4 the inhibition group; the continuous version used
  iCMS2-up scores as the state variable.
- **Robustness.** Seven perturbation scenarios — excluding three of the nine
  training datasets, removing 30% of the genes, deleting 30% of the enrichment
  methods, dropping four of the twelve combination algorithms, and adding
  Gaussian noise (0–100%, additive and multiplicative, extreme settings
  included) — were scored by the Pearson correlation between the perturbed and
  the original `DI`. MetaScreener outperformed conventional single-dataset,
  single-enrichment applications, including under extreme noise.
- **Independent validation.** The prioritised directional signatures were
  confirmed in 22 independent CRC datasets and in a massive-scale single-cell
  RNA-seq dataset.
- **Application.** In a signature-guided drug screen, parthenolide was
  identified and experimentally validated as a dual inhibitor of Wnt signalling
  and CRC metastasis.

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
install.packages("remotes")

remotes::install_github("IOBR/IOBR")          # not on CRAN or Bioconductor
BiocManager::install(c("GSVA", "GSEABase", "AUCell", "singscore", "UCell",
                       "decoupleR", "clusterProfiler", "limma", "fgsea"),
                     update = FALSE)
install.packages(c("metap", "doBy", "dplyr", "plyr", "reshape2"))

remotes::install_github("FangZY-Lab/MetaScreener")
library(MetaScreener)
```

The first installation takes 15–40 minutes because of the Bioconductor
dependency stack. New to the package? Start with the
**[step-by-step tutorial](TUTORIAL.md)**.

## Quick start

A reduced run — two of the 23 enrichment methods — that finishes in about
twenty seconds and shows the shape of the output:

```r
library(MetaScreener)

data("example_geneset")      # 10 Wnt signatures
data("example_diffgroup")    # per-dataset comparison table
data(list = c("diffGSE33113", "diffGSE33113_G",
              "diffGSE35896", "diffGSE35896_G",
              "diffGSE13067", "diffGSE13067_G"))

res <- DiffMetaScreener(
  expression_accession_vector = c("diffGSE33113", "diffGSE35896", "diffGSE13067"),
  comparisons                 = example_diffgroup,
  geneSets_gmt                = example_geneset,
  enrichment_method           = c("zscore_t", "zscore_wilcoxon"),
  min.sz = 2,
  max.sz = 10000
)

res$DI <- res$ADI + res$IDI
res[order(-res$DI), ]
```

```text
   signatures               ADI            IDI             DI
1  FLIER_WNT        18.74944260  -6.136770e-06    18.74943646
2  BIOCARTA_WNT     12.01998509  -2.179156e-02    11.99819353
3  HALLMARK_WNT     11.14547430  -2.411874e-03    11.14306243
4  WP_WNT            7.31220985  -3.819236e-02     7.27401749
5  KEGG_WNT          6.16164280  -2.953711e-01     5.86627170
6  FEVR_WNT_UP       5.30784588  -3.516706e-02     5.27267882
7  REACTOME_WNT      1.23345962  -5.840720e-02     1.17505242
8  GOBP_WNT          0.57010137  -1.005907e+00    -0.43580563
9  PID_WNT           0.09815316  -1.713012e+00    -1.61485884
10 FEVR_WNT_DN       0.01452660  -8.911991e+00    -8.89746440
```

`FLIER_WNT` carries activation evidence across all three cohorts, `FEVR_WNT_DN`
is a pure inhibition signature, and `GOBP_WNT` / `PID_WNT` carry both activation
and inhibition evidence with inhibition dominating, which is what the negative
`DI` reports. Run the full 23-method set for a result comparable to the
publication.

## Input contract

Objects are looked up **by name** in the calling environment:

| Object | Content |
| --- | --- |
| `<id>` | omics matrix, features × samples, calibrated (e.g. log2), unique feature IDs |
| `<id>_G` (Diff) | data frame with `Tag` and `group` |
| `<id>_V` (Cor) | data frame with `Tag` and `value` |

`Tag` must match `colnames(<id>)` — same values in the same order. Gene sets are
supplied as a `GeneSetCollection`, for example `GSEABase::getGmt("sets.gmt")`.

## A note on small k

`meanp` and `meanz` require at least four valid p-values. They are applied inside
each dataset; if you restrict `p_combine_method` to them while using few methods
or few datasets, they return `NA`, the pipeline substitutes 1, and every ADI/IDI
collapses to zero. The functions warn when this is about to happen; keep the
default (all twelve algorithms) to avoid it.

## Bundled example data

| Object | Dimensions | Content |
| --- | --- | --- |
| `diffGSE33113`, `diffGSE35896`, `diffGSE13067` | 21,502 genes × 80 / 54 / 58 samples | log2 expression, colorectal cancer |
| `<id>_G` | samples × 2 | `Tag`, `group` (CMS2 vs Other) |
| `corGSE33113`, `corGSE35896`, `corGSE13067` | 21,502 genes × 90 / 62 / 74 samples | log2 expression, colorectal cancer |
| `<id>_V` | samples × 2 | `Tag`, `value` |
| `example_diffgroup` | 3 × 5 | comparison table |
| `example_geneset` | `GeneSetCollection` of 10 | Wnt signatures |

The three cohorts are public GEO series (GSE13067, GSE33113, GSE35896). See
`?MetaScreener-data` for details.

## Citation

If MetaScreener contributes to your work, please cite:

> Zhao D, Zhao G, Yao M, Wu J, Fang Z. MetaScreener: a robust dual-mode framework
> for directional prioritization of actionable signatures through multi-dataset
> and multi-approach integration. *J Transl Med*. 2026;24(1):583.
> doi:10.1186/s12967-026-08019-y. PMID: 41845481.

```bibtex
@article{zhao2026metascreener,
  title   = {MetaScreener: a robust dual-mode framework for directional
             prioritization of actionable signatures through multi-dataset
             and multi-approach integration},
  author  = {Zhao, D. and Zhao, G. and Yao, M. and Wu, J. and Fang, Z.},
  journal = {Journal of Translational Medicine},
  year    = {2026},
  volume  = {24},
  number  = {1},
  pages   = {583},
  doi     = {10.1186/s12967-026-08019-y}
}
```

## Documentation

- `?DiffMetaScreener`, `?CorMetaScreener` — argument reference and worked examples
- `?MetaScreener-data` — the bundled cohorts
- [TUTORIAL.md](TUTORIAL.md) — a ten-minute walkthrough for first-time users
- [NEWS.md](NEWS.md) — release notes

## Contact

Dingkang Zhao — <dingkang.25@intl.zju.edu.cn>

Bug reports and feature requests:
<https://github.com/FangZY-Lab/MetaScreener/issues>.

## License

MIT © Dingkang Zhao. See [LICENSE](LICENSE).
