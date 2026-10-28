import pandas as pd
import pytest
from ucsd_pfs.dataset import COLUMN_MAP, audit_version_permutation, build_interim, clean_clinical


def _raw(n=3, **overrides):
    raw = {c: ["x"] * n for c in COLUMN_MAP}
    raw["ID"] = [f"UCSD-PTGBM-000{i + 1}_01" for i in range(n)]
    for c in COLUMN_MAP:
        if c.startswith("Days") or c in ("Overall survival", "Progression Free survival"):
            raw[c] = [str(10 * (i + 1)) for i in range(n)]
    raw["Patient's Age"] = [50, 60, 70][:n]
    raw.update(overrides)
    # Source sheet has trailing spaces on some headers.
    return pd.DataFrame(raw).rename(columns={"Surgery": "Surgery "}).rename(columns=str.strip)


def test_clean_clinical_renames_parses_and_coerces():
    raw = _raw(
        **{
            "Overall survival": ["Unknown", "12", "30"],
            "Avastin® (bevacizumab)": ["N", " No ", "Yes"],
        }
    )
    df = clean_clinical(raw)
    assert set(COLUMN_MAP.values()) <= set(df.columns)
    assert df["patient_id"].tolist() == ["UCSD-PTGBM-0001", "UCSD-PTGBM-0002", "UCSD-PTGBM-0003"]
    assert df["tp"].tolist() == [1, 1, 1]
    assert pd.isna(df.loc[0, "os_days_provided"]) and df.loc[1, "os_days_provided"] == 12
    assert df["bevacizumab"].tolist() == ["No", "No", "Yes"]


def test_clean_clinical_rejects_bad_ids():
    with pytest.raises(ValueError):
        clean_clinical(_raw(ID=["bad", "UCSD-PTGBM-0002_01", "UCSD-PTGBM-0003_01"]))


def test_audit_detects_row_permutation_and_interim_keeps_v1_survival():
    v1 = clean_clinical(_raw())
    # v3: every column but id/age moved one row down (what the v3 "ID correction" did).
    v3 = v1.copy()
    moved = [c for c in v1.columns if c not in ("id", "patient_id", "tp", "brats_id", "age")]
    for c in moved:
        v3[c] = v1[c].iloc[[1, 2, 0]].to_numpy()
    audit = audit_version_permutation(v1, v3).set_index("column")
    assert audit.loc["d_death", "moved_with_dates"] == 1.0
    assert audit.loc["age", "unmoved"] == 1.0

    neg = pd.DataFrame({"ID": ["UCSD-PTGBM-0002_01 "], "Dx": ["PsP"]})
    df = build_interim(v1, v3, neg)
    assert (df["os_days_provided"].to_numpy() == v1["os_days_provided"].to_numpy()).all()
    assert df["has_tumor"].tolist() == [True, False, True]
