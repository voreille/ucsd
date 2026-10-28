"""Build cleaned OS/PFS targets for the survival subset.

Always take survival targets from ``data/processed/survival_patients.parquet`` (built here), never
from raw spreadsheet columns. See references/survival_cleaning.md for the rationale.
"""

from typing import Literal

import numpy as np
import pandas as pd
import typer
from loguru import logger

from ucsd_pfs.config import (
    INTERIM_TIMEPOINTS,
    OS_MATCH_TOL_DAYS,
    POSTOP_WINDOW_DAYS,
    PROCESSED_SURVIVAL,
    PROCESSED_TIMEPOINTS,
)
from ucsd_pfs.timeline import (
    acquisition_order_consistent,
    add_timepoint_times,
    first_scan_in_window,
    patient_timeline,
)

app = typer.Typer()

PfsMode = Literal["events_only", "censor_at_last_followup"]


def recompute_os(pt: pd.DataFrame) -> pd.DataFrame:
    """OS from dates: initial surgery -> death (event) or -> last follow-up (censored)."""
    pt = pt.copy()
    anchored = pt["has_surgery_anchor"]
    pt["os_days_recomputed"] = np.where(pt["death_recorded"], pt["t_death"], pt["t_last_followup"])
    pt.loc[~anchored, "os_days_recomputed"] = np.nan
    pt["os_event"] = pt["death_recorded"].astype(int)
    return pt


def select_paper_subset(tp: pd.DataFrame) -> pd.DataFrame:
    """Rows the authors attached OS/PFS to, one per patient.

    Prefer the earliest such scan inside the post-op window, else the earliest one.
    """
    sub = tp[tp["os_days_provided"].notna() | tp["pfs_days_provided"].notna()].copy()
    sub["_outside"] = ~sub["days_since_surgery"].between(0, POSTOP_WINDOW_DAYS)
    sub = sub.sort_values(
        ["patient_id", "_outside", "acq_rel_last_followup", "tp"], na_position="last"
    ).drop(columns="_outside")
    n_conflict = (
        sub.groupby("patient_id")[["os_days_provided", "pfs_days_provided"]].nunique() > 1
    ).any(axis=1)
    if n_conflict.any():
        logger.warning(
            f"Patients with conflicting OS/PFS across rows: {list(n_conflict[n_conflict].index)}"
        )
    return sub.groupby("patient_id", as_index=False).head(1)


def qc_flags(row: pd.Series, tol: int) -> str:
    flags = []
    if not row["timeline_consistent"]:
        flags.append("timeline_inconsistent")
    if pd.isna(row["os_days_provided"]):
        flags.append("os_missing")
    if pd.isna(row["pfs_days_provided"]):
        flags.append("pfs_missing")
    if row["pfs_days_provided"] < 0:
        flags.append("pfs_negative")
    if row["pfs_days_provided"] > row["os_days_provided"]:
        flags.append("pfs_gt_os")
    if row["pfs_days_provided"] == row["os_days_provided"]:
        flags.append("pfs_eq_os")
    if pd.notna(row["os_days_recomputed"]) and pd.notna(row["os_days_provided"]):
        if abs(row["os_days_recomputed"] - row["os_days_provided"]) > tol:
            flags.append("os_mismatch_dates")
    if not row["has_surgery_anchor"]:
        flags.append("no_surgery_date")
    if not row["death_recorded"]:
        flags.append("no_death_date")
    if row["days_since_surgery"] > POSTOP_WINDOW_DAYS:
        flags.append("scan_outside_window")
    return ";".join(flags)


def select_cohort(tp: pd.DataFrame) -> pd.DataFrame:
    """Survival cohort: the provided-OS/PFS rows, plus patients the paper's rule adds.

    Paper rule: first scan within POSTOP_WINDOW_DAYS after initial surgery. It reproduces every
    patient of the provided subset and adds a few for whom only dates (no OS/PFS) are given.
    """
    sub = select_paper_subset(tp).assign(in_provided_subset=True)
    rule = first_scan_in_window(tp, POSTOP_WINDOW_DAYS, "days_since_surgery")
    extra = rule[~rule["patient_id"].isin(sub["patient_id"])].assign(in_provided_subset=False)
    return pd.concat([sub, extra], ignore_index=True).sort_values("patient_id")


def build_survival_table(tp: pd.DataFrame, tol: int = OS_MATCH_TOL_DAYS) -> pd.DataFrame:
    pt = recompute_os(patient_timeline(tp, tol))
    keep = [
        "patient_id",
        "id",
        "tp",
        "in_provided_subset",
        "t_acq",
        "t_acq_implied_by_os",
        "days_since_surgery",
        "days_since_last_surgery",
        "os_days_provided",
        "pfs_days_provided",
        "v3_os_days_provided",
        "v3_pfs_days_provided",
    ]
    df = select_cohort(tp)[keep].rename(columns={"id": "cohort_tp_id", "tp": "cohort_tp"})
    df = df.merge(pt, on="patient_id", how="left", validate="1:1")
    # Provided (v1) OS where available: it agrees with the dates and also covers undated surgery.
    provided = df["os_days_provided"].notna()
    df["os_days"] = df["os_days_provided"].fillna(df["os_days_recomputed"])
    df["os_days_source"] = np.where(
        provided, "provided_v1", np.where(df["os_days"].notna(), "recomputed_dates", "none")
    )
    df["pfs_days"] = df["pfs_days_provided"]
    df["qc_flags"] = df.apply(qc_flags, axis=1, tol=tol)
    return df.reset_index(drop=True)


def pfs_target(df: pd.DataFrame, mode: PfsMode) -> pd.DataFrame:
    """(time, event) for PFS under an explicit assumption; the data has no PFS event flag.

    events_only: every provided PFS is an observed progression; rows without PFS are dropped.
    censor_at_last_followup: rows without PFS are censored at last follow-up (from surgery).
    """
    out = df[["patient_id"]].copy()
    if mode == "events_only":
        out["pfs_time"] = df["pfs_days"]
        out["pfs_event"] = 1
        return out[df["pfs_days"].notna()].reset_index(drop=True)
    if mode == "censor_at_last_followup":
        known = df["pfs_days"].notna()
        out["pfs_time"] = df["pfs_days"].where(known, df["t_last_followup"])
        out["pfs_event"] = known.astype(int)
        return out[out["pfs_time"].notna()].reset_index(drop=True)
    raise ValueError(f"Unknown PFS mode: {mode}")


def summarize(tp: pd.DataFrame, surv: pd.DataFrame, tol: int) -> None:
    order_ok = acquisition_order_consistent(tp)
    logger.info(
        f"Timepoint suffix follows acquisition order for {order_ok.sum()}/{len(order_ok)} patients"
    )
    both = surv["os_days_recomputed"].notna() & surv["os_days_provided"].notna()
    dead = both & surv["death_recorded"]
    agree = (surv["os_days_recomputed"] - surv["os_days_provided"]).abs() <= tol
    logger.info(
        f"Survival cohort: {len(surv)} patients (paper: 94), of which "
        f"{surv['in_provided_subset'].sum()} have provided OS/PFS"
    )
    logger.info("OS source: " + str(surv["os_days_source"].value_counts().to_dict()))
    logger.info(
        f"Provided OS == recomputed OS (±{tol} d): {(agree & dead).sum()}/{dead.sum()} "
        "with death date"
    )
    # The same check with the unrepaired v3 OS column shows why v1 OS is used.
    v3_agree = (surv["os_days_recomputed"] - surv["v3_os_days_provided"]).abs() <= tol
    logger.info(f"v3 OS == recomputed OS (±{tol} d): {(v3_agree & surv['death_recorded']).sum()}")
    delay = surv["days_since_surgery"].describe().round(1).to_dict()
    logger.info(f"Scan delay after surgery (days): {delay}")
    in_sub = set(surv.loc[surv["in_provided_subset"], "patient_id"])
    for since in ("days_since_surgery", "days_since_last_surgery"):
        in_rule = set(first_scan_in_window(tp, POSTOP_WINDOW_DAYS, since)["patient_id"])
        logger.info(
            f"Rule cohort (first scan with {since} in [0, {POSTOP_WINDOW_DAYS}]): {len(in_rule)} "
            f"patients; overlap with provided subset {len(in_rule & in_sub)}; "
            f"union {len(in_rule | in_sub)}"
        )
    flags = surv["qc_flags"].str.split(";").explode()
    logger.info("QC flags:\n" + flags[flags != ""].value_counts().to_string())


@app.command()
def main(tol: int = OS_MATCH_TOL_DAYS):
    tp = add_timepoint_times(pd.read_parquet(INTERIM_TIMEPOINTS))
    surv = build_survival_table(tp, tol)
    summarize(tp, surv, tol)
    PROCESSED_SURVIVAL.parent.mkdir(parents=True, exist_ok=True)
    tp.to_parquet(PROCESSED_TIMEPOINTS, index=False)
    surv.to_parquet(PROCESSED_SURVIVAL, index=False)
    logger.success(f"Wrote {PROCESSED_TIMEPOINTS} and {PROCESSED_SURVIVAL}")


if __name__ == "__main__":
    app()
