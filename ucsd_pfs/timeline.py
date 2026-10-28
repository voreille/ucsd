"""Convert per-acquisition day offsets into per-patient timelines.

All ``d_*`` columns are days relative to the acquisition date of that row's MRI. Here they are
re-expressed as days since the patient's initial surgery/treatment/diagnosis (``t_*`` columns,
day 0 = initial surgery), which is the origin of OS and PFS in the data dictionary.
"""

import numpy as np
import pandas as pd

from ucsd_pfs.config import OS_MATCH_TOL_DAYS


def add_timepoint_times(df: pd.DataFrame) -> pd.DataFrame:
    """Add t_* columns (days since initial surgery) for each timepoint row."""
    df = df.copy()
    df["t_acq"] = -df["d_initial_surgery"]
    for src, dst in [
        ("d_last_surgery", "t_last_surgery"),
        ("d_death", "t_death"),
        ("d_last_followup", "t_last_followup"),
        ("d_first_radiation", "t_first_radiation"),
        ("d_first_chemo", "t_first_chemo"),
    ]:
        df[dst] = df["t_acq"] + df[src]
    # Days between the scan and the most recent surgery before it.
    df["days_since_last_surgery"] = -df[["d_initial_surgery", "d_last_surgery"]].max(axis=1)
    # Scan time implied by the provided OS (OS is measured from initial surgery to death).
    df["t_acq_implied_by_os"] = df["os_days_provided"] - df["d_death"]
    # Scan delay after initial surgery, falling back on the OS-implied value when undated.
    df["days_since_surgery"] = df["t_acq"].fillna(df["t_acq_implied_by_os"])
    # Acquisition order that does not need a surgery date (last follow-up is a per-patient date).
    df["acq_rel_last_followup"] = -df["d_last_followup"]
    return df


def _spread(s: pd.Series) -> float:
    s = s.dropna()
    return float(s.max() - s.min()) if len(s) > 1 else 0.0


def patient_timeline(tp: pd.DataFrame, tol: int = OS_MATCH_TOL_DAYS) -> pd.DataFrame:
    """One row per patient: anchored death/follow-up dates and cross-timepoint consistency."""
    g = tp.groupby("patient_id")
    out = pd.DataFrame(
        {
            "n_timepoints": g.size(),
            "has_surgery_anchor": g["t_acq"].apply(lambda s: s.notna().any()),
            "t_death": g["t_death"].median(),
            "t_last_followup": g["t_last_followup"].max(),
            "death_recorded": g["d_death"].apply(lambda s: s.notna().any()),
            # Same calendar date seen from different acquisitions should agree.
            "anchor_spread_days": g.apply(
                lambda d: _spread(d["d_initial_surgery"] - d["d_last_followup"]),
                include_groups=False,
            ),
            "death_spread_days": g["t_death"].apply(_spread),
        }
    )
    out["timeline_consistent"] = (out["anchor_spread_days"] <= tol) & (
        out["death_spread_days"] <= tol
    )
    return out.reset_index()


def first_scan_in_window(
    tp: pd.DataFrame, window_days: int, since: str = "days_since_last_surgery"
) -> pd.DataFrame:
    """Per patient, the earliest timepoint acquired within [0, window_days] of surgery."""
    cand = tp[tp[since].between(0, window_days)]
    cand = cand.sort_values(["patient_id", "acq_rel_last_followup", "tp"], na_position="last")
    return cand.groupby("patient_id", as_index=False).head(1)


def acquisition_order_consistent(tp: pd.DataFrame) -> pd.Series:
    """Per patient: does the _NN timepoint suffix follow acquisition order?"""

    def ok(d: pd.DataFrame) -> bool:
        a = d.sort_values("tp")["acq_rel_last_followup"].dropna().to_numpy()
        return bool(np.all(np.diff(a) >= 0))

    return tp.groupby("patient_id").apply(ok, include_groups=False)
