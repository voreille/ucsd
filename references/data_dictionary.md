# Clinical data dictionary (cleaned)

Source sheet: `Cinical info` (sic) in the TCIA clinical xlsx. There is one row per MRI timepoint, with IDs like `UCSD-PTGBM-XXXX_NN`. Cleaning is done in `ucsd_pfs/dataset.py`.

Raw quirks that cleaning handles:
- Some headers end in a space (`"Surgery "`, `"Days from Acquisition to Date of initial surgery, treatment or diagnosis "`, ...).
- The `Data Dictionary` sheet is padded with about 1M empty rows. Don't read it whole.
- `Unknown` appears inside numeric columns.
- Typos and duplicate categories: `Gliobastoma`, `"Asian "`, bevacizumab `N` vs `No`.
- Typos in headers: "Acquistion", "Surgery extend".

| cleaned name | raw header | notes |
|---|---|---|
| id, patient_id, tp | ID | `patient_id` = `UCSD-PTGBM-XXXX`, `tp` = NN; tp order follows acquisition order (177/178) |
| brats_id | BraTS ID | BraTS 2024 subject ID |
| field_strength, manufacturer, scanner_model | Magnetic Field Strength, Manufacturer, Manufacturer's Model Name | all 3T GE; MR750 or Signa HDxt |
| age | Patient's Age | at scan, in years |
| sex, race, ethnicity | Sex at birth, Race, Ethnicity | |
| primary_diagnosis, who2021_diagnosis, non_who2021_diagnosis, grade | | |
| mgmt, idh, codel_1p19q, atrx | MGMT, IDH, 1p19q, ATRX | lots missing |
| d_initial_surgery | Days from Acquisition to Date of initial surgery, treatment or diagnosis | anchor for OS/PFS; missing when surgery = No |
| d_last_surgery | ... last surgery prior to scan if different | |
| d_last_followup, d_death | ... last follow-up / death | |
| os_days_provided, pfs_days_provided | Overall survival, Progression Free survival | **taken from v1** (see survival_cleaning.md) |
| v3_os_days_provided, v3_pfs_days_provided | same, as published in v3 | mis-assigned; audit only |
| surgery, n_surgeries, extent_of_resection | Surgery, Number of surgeries, Surgery extend | GTR / STR / Biopsy |
| radiation, n_radiation_courses, d_first_radiation, d_last_radiation | | |
| first_chemo, d_first_chemo | 1st Chemo type, ... | Temodar / CCNU |
| bevacizumab, d_last_bevacizumab | Avastin® (bevacizumab), ... | |
| other_treatment, d_other_treatment | | dates are always Unknown |
| negative_case_dx, has_tumor | `details_of_negative_cases_TCIA.xlsx` | 51 timepoints that only show post-treatment change |

Derived columns (`timeline.py`), all in days since initial surgery:
- `t_acq`, `t_death`, `t_last_followup`, `t_last_surgery`, `t_first_radiation`, `t_first_chemo`.
- `days_since_surgery`: uses `t_acq`, falling back on `os_days_provided − d_death`.
- `days_since_last_surgery`.
- `acq_rel_last_followup`: for ordering scans when there's no surgery date.
