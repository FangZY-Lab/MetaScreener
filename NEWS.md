# MetaScreener 1.1.0

## New documentation

* `man/figures/MetaScreener-workflow.svg` — a schematic of the eight-stage
  pipeline, referenced from the README.
* `TUTORIAL.md` — a step-by-step walkthrough for first-time users: install,
  load the bundled data, run a first screen, interpret ADI/IDI, adapt the code
  to your own cohorts, choose methods, and troubleshoot.
* The README was rewritten around the statistics of the method (why the two
  directions are kept apart, what the indices mean, the method space) rather
  than around a list of function calls.

## Bug fixes

* The two exported functions call `library()` on 24 packages but only 15 were
  declared in `DESCRIPTION`. On a clean machine `caret`, `coin`,
  `correlation`, `GSA`, `stringr`, `tidyverse`, `viper` and `wdm` were missing
  and the first call failed with "there is no package called ...". They are now
  declared, so `install_github()` pulls them in.
* `library(devtools)` was called although no devtools function is used; the
  call was removed and the package no longer needs it at runtime.
* Both functions now validate their inputs before scoring: a missing dataset
  object, a missing `<id>_G` / `<id>_V` table, a wrong column name, or a `Tag`
  column that does not match `colnames()` produces a message naming the
  offending objects instead of an obscure failure deep inside the pipeline.
* `meanp` and `meanz` require at least four valid p-values at each of the two
  combination steps. Below that they return `NA`, which the pipeline replaced
  with 1, silently flattening every ADI and IDI to zero. The functions now warn
  before this happens.

## Documentation of the interface

* `DiffMetaScreener()` and `CorMetaScreener()` gained `\value{}` and runnable
  `\examples{}` sections (`man/*.Rd` regenerated with roxygen2), so `?` now
  documents the returned data frame and shows a working call.

## Metadata

* Contact address updated to <dingkang.25@intl.zju.edu.cn>.
* `Title` and `Description` rewritten to state the statistical contribution.

No analysis method was changed.

---

# MetaScreener 1.0.1

## Packaging

* The example datasets are now ordinary files in `data/` instead of Git LFS
  pointers, so they survive `remotes::install_github("FangZY-Lab/MetaScreener")`.
  Previously the installed package contained fourteen 128-byte pointer files.
* Renamed the data files from `.Rdata` to `.rda`. `data()` only recognises the
  extensions `.rda`, `.RData` and `.rdata`, so the old names were invisible to
  it even when the content was correct; `data(diffGSE33113)` now works.
* `LazyData: true`, so the bundled objects can also be reached directly, e.g.
  `MetaScreener::example_geneset`.
* `DESCRIPTION` declares the MIT licence (`License: MIT + file LICENSE`) and
  adds `URL` / `BugReports`. The licence placeholder left by `usethis` is gone.
* Added `LICENSE`, `LICENSE.md`, `NEWS.md` and `?MetaScreener-data`, which
  documents all fourteen bundled objects.
* Removed `.Rhistory` and `data/desktop.ini` from the repository, dropped the
  LFS rules in `.gitattributes`, and expanded `.gitignore` and `.Rbuildignore`.
* README rewritten around the `data()` workflow.

## Documentation

* The README now records that `meanp` and `meanz` need at least four valid
  p-values at both combining steps. Passing one of them with few methods and
  few datasets silently flattens ADI and IDI to zero.

No analysis code changed.
