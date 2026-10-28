import numpy as np
import pandas as pd
from ucsd_pfs.survival import build_survival_table, pfs_target
from ucsd_pfs.timeline import add_timepoint_times, first_scan_in_window, patient_timeline


def _tp():
    # Patient A: surgery day 0, scans at day 30 and 400, dies day 600 (provided OS on scan 1).
    # Patient B: surgery day 0, scan at day 20, alive at day 900, no provided OS/PFS.
    rows = [
        dict(
            id="P-A_01",
            patient_id="A",
            tp=1,
            d_initial_surgery=-30,
            d_last_surgery=np.nan,
            d_death=570,
            d_last_followup=560,
            os_days_provided=600,
            pfs_days_provided=200,
        ),
        dict(
            id="P-A_02",
            patient_id="A",
            tp=2,
            d_initial_surgery=-400,
            d_last_surgery=-100,
            d_death=200,
            d_last_followup=190,
            os_days_provided=np.nan,
            pfs_days_provided=np.nan,
        ),
        dict(
            id="P-B_01",
            patient_id="B",
            tp=1,
            d_initial_surgery=-20,
            d_last_surgery=np.nan,
            d_death=np.nan,
            d_last_followup=880,
            os_days_provided=np.nan,
            pfs_days_provided=np.nan,
        ),
    ]
    df = pd.DataFrame(rows)
    for c in ("d_first_radiation", "d_first_chemo", "v3_os_days_provided", "v3_pfs_days_provided"):
        df[c] = np.nan
    return add_timepoint_times(df)


def test_timeline_anchors_on_initial_surgery():
    tp = _tp()
    assert tp["t_acq"].tolist() == [30, 400, 20]
    assert tp["t_death"].tolist()[:2] == [600, 600]
    assert tp["days_since_last_surgery"].tolist() == [30, 100, 20]
    pt = patient_timeline(tp).set_index("patient_id")
    assert pt["timeline_consistent"].all()
    assert pt.loc["A", "t_last_followup"] == 590


def test_first_scan_in_window_one_per_patient():
    sel = first_scan_in_window(_tp(), 90)
    assert sel["id"].tolist() == ["P-A_01", "P-B_01"]
    assert sel["days_since_last_surgery"].between(0, 90).all()


def test_survival_table_os_and_censoring():
    s = build_survival_table(_tp()).set_index("patient_id")
    assert s.loc["A", "os_days"] == 600 and s.loc["A", "os_event"] == 1
    assert s.loc["A", "os_days_source"] == "provided_v1" and s.loc["A", "in_provided_subset"]
    # B: no provided OS, alive -> censored at last follow-up, recomputed from dates.
    assert s.loc["B", "os_days"] == 900 and s.loc["B", "os_event"] == 0
    assert s.loc["B", "os_days_source"] == "recomputed_dates"
    assert not s[["os_days", "os_event"]].isna().any().any()


def test_pfs_target_modes():
    s = build_survival_table(_tp())
    ev = pfs_target(s, "events_only")
    assert ev["patient_id"].tolist() == ["A"] and ev["pfs_event"].tolist() == [1]
    cens = pfs_target(s, "censor_at_last_followup").set_index("patient_id")
    assert cens.loc["B", "pfs_time"] == 900 and cens.loc["B", "pfs_event"] == 0
