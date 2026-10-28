# OS/PFS cleaning log

Sources:
- clinical spreadsheet v1 (`UCSD_PTGBM-clinical-information_v1_2025-25-03.xlsx`)
- clinical spreadsheet v3 (`..._v3_2026-12-Mar.xlsx`)
- Gagnon et al., *Sci Data* 2026, doi 10.1038/s41597-025-06499-z

Code: `ucsd_pfs/dataset.py`, `timeline.py`, `survival.py`, `reproduce.py`. Walkthrough: `notebooks/0.01-vo-clinical-audit.ipynb`.

## Definitions (data dictionary)
- All `Days from Acquisition to ...` columns are day offsets from **that row's MRI acquisition date**. Negative means before the scan.
- OS = days from initial surgery to death. PFS = days from initial surgery to first sign of progression.
- There is **no event or censoring column** for either OS or PFS. A missing value is the string `Unknown`.

## Finding 1: v3 is a row permutation of v1, and OS/PFS should not have moved
- TCIA says v3 "corrected patient IDs incorrectly assigned" in earlier versions.
- Each of the 243 v3 rows matches exactly one v1 row on all date offsets (2 matches are ambiguous).
- Every column moved with the dates, except `ID`, `BraTS ID` and `Patient's Age`, which stayed in place. Table: `reports/tables/v1_v3_column_audit.csv`.
- Consistency check, OS ≟ `d_death − d_initial_surgery` on the same ID (±3 days):

  | OS from | dates from | agree |
  |---|---|---|
  | v1 | v1 | 0/61 |
  | v3 | v3 | 0/61 |
  | **v1** | **v3** | **60/62** |

- The v3 dates hold together within each patient: the initial-surgery date seen from different scans agrees for 46/49 multi-scan patients, against 19/49 in v1.

**Conclusion.** In v1, OS/PFS were already attached to the correct ID and the other clinical columns were not. The v3 fix moved the clinical columns to the right rows but also moved OS/PFS away from them. **We use v3 for everything except OS/PFS, which come from v1 by ID.** The v3 values are kept as `v3_os_days_provided` / `v3_pfs_days_provided`, for audit only.

Corroboration:
- With the repaired data, every row carrying OS/PFS is a scan 0–90 days after initial surgery (median 28 days). That is exactly the paper's rule for the survival subset. With v3 as published, the median is 181 days and the maximum 2657.
- MGMT in the survival subset: the paper reports 17 methylated / 44 unmethylated. The repaired data gives 19/41; unrepaired v3 gives 22/19.

## Finding 2: survival cohort
- The provided subset (rows with v1 OS or PFS) has **85 patients**; the paper reports 94.
- Applying the paper's rule (first scan within 0–90 days of *initial* surgery) gives **96 patients**. All 85 provided patients are among them, plus 11 with dates but no OS/PFS. Measuring from the *last* surgery gives 90 patients, of which only 65 overlap, so the rule uses initial surgery.
- `survival_patients.parquet` holds all 96, with `in_provided_subset` marking the 85.
- Patients 0066 and 0101/0103/0109/0116 have OS/PFS on two rows (same values). We keep the earliest in-window row.

## Finding 3: OS
- `os_days` = provided v1 OS (81 patients). Otherwise it is recomputed from dates (14): death gives an event, else the patient is censored at the latest follow-up. `os_days_source` records which.
- `os_event` = 1 if a death date exists. 2 provided-OS patients have no death date (`no_death_date`).
- Almost every provided OS is a death. The authors' subset looks like deceased patients only, so censoring hardly appears in it.
- Two provided OS values disagree with the dates: 0075 (246 vs 231) and 0107 (364 vs 399). We keep the provided values and flag them `os_mismatch_dates`.
- 20 patients have no surgery date (`surgery = No`; the anchor is "initial treatment or diagnosis"). For them, the scan delay is taken from `OS − d_death`.
- 0004 and 0145 have inconsistent anchors across their scans (`timeline_inconsistent`). For 0145 the gap is 3652 days, which looks like a 10-year typo.

## Finding 4: PFS (unresolved)
- PFS **cannot be recomputed**: no progression date is given.
- PFS ≤ OS in every row.
- 0061 has PFS = −8 (`pfs_negative`).
- 3 patients have PFS = OS (`pfs_eq_os`): either death without documented progression, or a censored value.
- **Event status is unknown.** `survival.pfs_target(df, mode)` makes the assumption explicit:
  - `events_only`: every provided PFS counts as a progression; patients without PFS are dropped.
  - `censor_at_last_followup`: patients without PFS are censored at their last follow-up.

## Finding 5: reproduction of the paper (`make reproduce`)
- **Table 1, total column:**
  - Matches: 243 timepoints, 178 patients, age 55.7 ± 13.4, MGMT 76/73.
  - Doesn't match: sex is 129 M / 49 F (paper 116/62; its own text says 121 men), and GTR is 83 (paper 87). v1 gives 119/59 and 90, which doesn't match either. Unresolved.
- **Fig. 2** (paper: log-rank p = 0.00002; "survival" taken to mean OS):

  | variant | n | p |
  |---|---|---|
  | provided subset, v1 OS (all treated as events) | 56 | 0.32 |
  | unrepaired v3 OS | 39 | 0.36 |
  | rule cohort (96), cleaned OS | 66 | 5e-4 |
  | all patients, OS recomputed from dates | 91 | 3e-7 |

  None reproduces p = 2e-5 exactly. Using the whole cohort with OS from dates gives the expected direction and significance.

## Open questions (maybe ask the authors / TCIA)
1. Confirm that the v3 OS/PFS columns were moved by mistake.
2. Which 94 patients, which OS definition and which censoring were used in Table 1 / Fig. 2?
3. Is a missing PFS a censored patient or missing data? How should PFS = OS be read?
4. Which sex/EOR values are correct, given the Table 1 discrepancy?
