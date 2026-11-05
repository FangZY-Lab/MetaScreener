# MetaScreener 1.0.1

## Packaging

* The example datasets are now ordinary files in `data/` instead of Git LFS
  pointers, so `data(diffGSE33113)` works after
  `remotes::install_github("FangZY-Lab/MetaScreener")`. Previously the
  installed package contained 128-byte pointer files and every `data()` call
  failed.
* `LazyData: true`, so the bundled objects can also be referenced directly.
* `DESCRIPTION` declares the MIT licence (`License: MIT + file LICENSE`) and
  adds `URL` / `BugReports`. The licence placeholder left by `usethis` is gone.
* Added `LICENSE`, `LICENSE.md` and `?MetaScreener-data`, which documents all
  fourteen bundled objects.
* Removed `.Rhistory` and `data/desktop.ini` from the repository, dropped the
  LFS rules in `.gitattributes`, and expanded `.gitignore` and `.Rbuildignore`.
* README rewritten around the `data()` workflow.

No analysis code changed.
