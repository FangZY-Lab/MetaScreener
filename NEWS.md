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
