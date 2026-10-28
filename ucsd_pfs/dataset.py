"""Download and clean the UCSD-PTGBM clinical spreadsheets.

The interim table is v3 for every column except the provided OS/PFS, which are taken from v1
(see references/survival_cleaning.md for why). The v3 OS/PFS values are kept, prefixed with
``v3_``, for auditing only.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import requests
import typer
from loguru import logger

from ucsd_pfs.config import (
    CLINICAL_FILES,
    CLINICAL_SHEET,
    INTERIM_TIMEPOINTS,
    NEGATIVE_CASES_FILE,
    RAW_DATA_DIR,
    TABLES_DIR,
    TCIA_UPLOADS,
)

app = typer.Typer()

# Raw column (after .strip()) -> snake_case name. "d_" columns are day offsets relative to the
# acquisition date of that row's MRI (negative = before the scan).
COLUMN_MAP = {
    "ID": "id",
    "BraTS ID": "brats_id",
    "Magnetic Field Strength": "field_strength",
    "Manufacturer": "manufacturer",
    "Manufacturer's Model Name": "scanner_model",
    "Patient's Age": "age",
    "Sex at birth": "sex",
    "Race": "race",
    "Ethnicity": "ethnicity",
    "Primary Diagnosis": "primary_diagnosis",
    "Days from Acquisition to Date of initial surgery, treatment or diagnosis": (
        "d_initial_surgery"
    ),
    "Days from Acquisition to Date of last surgery prior to scan if different": "d_last_surgery",
    "WHO 2021 Diagnosis": "who2021_diagnosis",
    "Non WHO 2021 Diagnosis": "non_who2021_diagnosis",
    "Grade": "grade",
    "MGMT": "mgmt",
    "IDH": "idh",
    "1p19q": "codel_1p19q",
    "ATRX": "atrx",
    "Days from Acquisition to Date of last follow-up": "d_last_followup",
    "Days from Acquisition to Date of death": "d_death",
    "Overall survival": "os_days_provided",
    "Progression Free survival": "pfs_days_provided",
    "Surgery": "surgery",
    "Number of surgeries": "n_surgeries",
    "Surgery extend": "extent_of_resection",
    "Radiation": "radiation",
    "Number of radiation courses": "n_radiation_courses",
    "Days from Acquistion to Date of first radiation": "d_first_radiation",
    "Days from Acquistion to Date of last radiation prior to scan": "d_last_radiation",
    "1st Chemo type": "first_chemo",
    "Days from Acquisition to Date of 1st chemo start": "d_first_chemo",
    "Avastin® (bevacizumab)": "bevacizumab",
    "Days from Acquisition to Date of last Avastin treatment prior to scan": "d_last_bevacizumab",
    "Other treatment prior to scan": "other_treatment",
    "Days from Acquisition to Other treatment dates": "d_other_treatment",
}

NUMERIC_COLUMNS = [
    "age",
    "grade",
    "n_surgeries",
    "n_radiation_courses",
    "os_days_provided",
    "pfs_days_provided",
] + [c for c in COLUMN_MAP.values() if c.startswith("d_")]

# Category fixes: strip whitespace first, then apply these.
CATEGORY_FIXES = {
    "who2021_diagnosis": {"Gliobastoma": "Glioblastoma"},
    "bevacizumab": {"N": "No"},
    "race": {"Native Hawaiian or Other Pacific Islander/ Asian": "Multiple"},
}

MISSING_TOKENS = {"Unknown", "unknown", "Unknown or not reported", "NA", "N/A", ""}

# Columns that v3 kept in place when re-assigning rows (see audit_version_permutation).
V1_COLUMNS_KEPT = ["os_days_provided", "pfs_days_provided"]


def download(filename: str, out_dir: Path = RAW_DATA_DIR, force: bool = False) -> Path:
    out = out_dir / filename
    if out.exists() and not force:
        logger.info(f"{out.name} already present")
        return out
    url = f"{TCIA_UPLOADS}/{filename}"
    logger.info(f"Downloading {url}")
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_bytes(r.content)
    return out


def read_clinical_raw(path: Path) -> pd.DataFrame:
    """Read the clinical sheet with stripped column names and raw values."""
    df = pd.read_excel(path, sheet_name=CLINICAL_SHEET)
    df.columns = [c.strip() for c in df.columns]
    return df


def clean_clinical(raw: pd.DataFrame) -> pd.DataFrame:
    """Rename to snake_case, coerce numerics, normalise categories, add patient/timepoint ids."""
    missing = set(COLUMN_MAP) - set(raw.columns)
    if missing:
        raise ValueError(f"Unexpected clinical sheet layout, missing columns: {missing}")
    df = raw.rename(columns=COLUMN_MAP)[list(COLUMN_MAP.values())].copy()

    for col in df.columns:
        if col in NUMERIC_COLUMNS:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        elif col not in ("id", "brats_id"):
            s = df[col].astype("string").str.strip()
            s = s.mask(s.isin(MISSING_TOKENS))
            df[col] = s.replace(CATEGORY_FIXES.get(col, {}))

    parts = df["id"].str.extract(r"^(?P<patient_id>UCSD-PTGBM-\d{4})_(?P<tp>\d{2})$")
    if parts.isna().any().any():
        raise ValueError("Unparseable IDs: " + ", ".join(df.loc[parts.isna().any(axis=1), "id"]))
    df.insert(1, "patient_id", parts["patient_id"])
    df.insert(2, "tp", parts["tp"].astype(int))
    return df.sort_values(["patient_id", "tp"]).reset_index(drop=True)


def audit_version_permutation(v1: pd.DataFrame, v3: pd.DataFrame) -> pd.DataFrame:
    """Describe how v3 rows relate to v1 rows, per column.

    Each v3 row is matched to the v1 row with identical date offsets. For every column we report
    the fraction of rows where v3 equals v1 at the same ID ("unmoved") and at the matched row
    ("moved_with_dates").
    """
    v1 = v1.set_index("id").loc[v3["id"]].reset_index()
    as_str = lambda d: d.astype(object).where(d.notna(), "NA").map(str)  # noqa: E731
    s1, s3 = as_str(v1), as_str(v3)
    dcols = [c for c in v3.columns if c.startswith("d_")]
    k1, k3 = s1[dcols].agg("|".join, axis=1), s3[dcols].agg("|".join, axis=1)
    lookup: dict[str, list[int]] = {}
    for i, k in k1.items():
        lookup.setdefault(k, []).append(i)
    match = np.array([lookup.get(k, [-1])[0] for k in k3])
    n_ambiguous = sum(len(lookup.get(k, [])) > 1 for k in k3)
    ok = match >= 0
    logger.info(
        f"v3 rows matched to a v1 row by dates: {ok.sum()}/{len(v3)} ({n_ambiguous} ambiguous)"
    )
    rows = []
    for c in v3.columns:
        rows.append(
            {
                "column": c,
                "unmoved": (s1[c].values == s3[c].values).mean(),
                "moved_with_dates": (s1[c].values[match[ok]] == s3[c].values[ok]).mean(),
            }
        )
    return pd.DataFrame(rows)


def build_interim(v1: pd.DataFrame, v3: pd.DataFrame, negatives: pd.DataFrame) -> pd.DataFrame:
    """v3 clinical data + v1 OS/PFS (by ID) + negative-case labels."""
    df = v3.rename(columns={c: f"v3_{c}" for c in V1_COLUMNS_KEPT})
    df = df.merge(v1[["id", *V1_COLUMNS_KEPT]], on="id", how="left", validate="1:1")
    neg = negatives.rename(columns={"ID": "id", "Dx": "negative_case_dx"})
    neg["id"] = neg["id"].str.strip()
    df = df.merge(neg[["id", "negative_case_dx"]], on="id", how="left", validate="1:1")
    df["has_tumor"] = df["negative_case_dx"].isna()
    return df


@app.command()
def main(force_download: bool = False):
    paths = {v: download(f, force=force_download) for v, f in CLINICAL_FILES.items()}
    neg_path = download(NEGATIVE_CASES_FILE, force=force_download)

    v1 = clean_clinical(read_clinical_raw(paths["v1"]))
    v3 = clean_clinical(read_clinical_raw(paths["v3"]))
    negatives = pd.read_excel(neg_path)

    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    audit = audit_version_permutation(v1, v3)
    audit.to_csv(TABLES_DIR / "v1_v3_column_audit.csv", index=False)
    logger.info("v1 -> v3 column audit:\n" + audit.round(2).to_string(index=False))

    df = build_interim(v1, v3, negatives)
    INTERIM_TIMEPOINTS.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(INTERIM_TIMEPOINTS, index=False)
    logger.success(
        f"Wrote {INTERIM_TIMEPOINTS} ({len(df)} timepoints, {df.patient_id.nunique()} patients)"
    )


if __name__ == "__main__":
    app()
