# ucsd-pfs

This project reproduces the UCSD-PTGBM dataset paper (Gagnon et al., *Sci Data* 2026, doi 10.1038/s41597-025-06499-z). After that, it predicts progression-free survival from the first post-surgery MRI timepoint, starting with clinical variables only.

Data: [UCSD-PTGBM on TCIA](https://www.cancerimagingarchive.net/collection/ucsd-ptgbm/) (doi 10.7937/fwv2-dt74, CC BY 4.0). Data is downloaded, never committed. The code is MIT.

> **Warning:** the OS/PFS columns of the published v3 clinical spreadsheet are attached to the wrong rows. This repo rebuilds them from v1. See [references/survival_cleaning.md](references/survival_cleaning.md).

## Quickstart
```bash
uv sync
make reproduce     # download -> clean -> data/processed/ -> reports/tables + reports/figures
make notebooks     # execute notebooks
make test lint
```

## Layout ([Cookiecutter Data Science v2](https://cookiecutter-data-science.drivendata.org/))
```
data/{raw,interim,processed,external}   downloaded / intermediate / final data (git-ignored)
notebooks/                              exploration, named <n.nn>-<initials>-<description>.ipynb
references/                             data dictionary, cleaning log
reports/{figures,tables}                generated outputs (git-ignored)
ucsd_pfs/                               source package
  config.py      paths and constants
  dataset.py     download + clean the clinical spreadsheets
  timeline.py    day offsets -> per-patient timeline anchored on initial surgery
  survival.py    survival cohort, cleaned OS/PFS targets, QC flags
  reproduce.py   paper Table 1 and Fig. 2
  plots.py       Kaplan-Meier helpers
  modeling/      PFS models (next step)
tests/
```
