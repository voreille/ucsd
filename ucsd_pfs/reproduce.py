"""Reproduce Table 1 and Fig. 2 of Gagnon et al., Sci Data 2026 (doi 10.1038/s41597-025-06499-z).

Every number is computed for several data variants and written next to the paper's value in
reports/tables/, so that the effect of the v3 OS/PFS mis-assignment and of each cohort choice is
explicit.
"""

import pandas as pd
import typer
from loguru import logger

from ucsd_pfs.config import FIGURES_DIR, PROCESSED_SURVIVAL, PROCESSED_TIMEPOINTS, TABLES_DIR
from ucsd_pfs.plots import km_by_group

app = typer.Typer()

PAPER_TABLE1 = {
    "total": {
        "timepoints": 243,
        "patients": 178,
        "age_mean": 56,
        "age_sd": 13,
        "male": 116,
        "female": 62,
        "gtr": 87,
        "mgmt_meth": 76,
        "mgmt_unmeth": 73,
    },
    "survival_subset": {
        "timepoints": 94,
        "patients": 94,
        "age_mean": 58,
        "age_sd": 13,
        "male": 67,
        "female": 27,
        "gtr": 38,
        "mgmt_meth": 17,
        "mgmt_unmeth": 44,
    },
}
PAPER_FIG2_P = 0.00002


def table1_column(tp: pd.DataFrame, per_patient: pd.DataFrame) -> dict:
    """Table 1 statistics. Counts per patient use ``per_patient`` (one row per patient)."""
    return {
        "timepoints": len(tp),
        "patients": per_patient["patient_id"].nunique(),
        "age_mean": round(per_patient["age"].mean(), 1),
        "age_sd": round(per_patient["age"].std(), 1),
        "male": int((per_patient["sex"] == "Male").sum()),
        "female": int((per_patient["sex"] == "Female").sum()),
        "gtr": int((per_patient["extent_of_resection"] == "GTR").sum()),
        "mgmt_meth": int((tp["mgmt"] == "Methylated").sum()),
        "mgmt_unmeth": int((tp["mgmt"] == "Unmethylated").sum()),
    }


def reproduce_table1(tp: pd.DataFrame, surv: pd.DataFrame) -> pd.DataFrame:
    first = tp.sort_values(["patient_id", "tp"]).groupby("patient_id", as_index=False).head(1)
    cohort_tp = tp[tp["id"].isin(surv["cohort_tp_id"])]
    provided_tp = tp[tp["id"].isin(surv.loc[surv["in_provided_subset"], "cohort_tp_id"])]
    v3_tp = tp[tp["v3_os_days_provided"].notna() | tp["v3_pfs_days_provided"].notna()]
    v3_tp = v3_tp.sort_values(["patient_id", "tp"]).groupby("patient_id", as_index=False).head(1)
    cols = {
        ("total", "paper"): PAPER_TABLE1["total"],
        ("total", "ours"): table1_column(tp, first),
        ("survival_subset", "paper"): PAPER_TABLE1["survival_subset"],
        ("survival_subset", "ours_provided_os"): table1_column(provided_tp, provided_tp),
        ("survival_subset", "ours_rule_cohort"): table1_column(cohort_tp, cohort_tp),
        ("survival_subset", "v3_unrepaired"): table1_column(v3_tp, v3_tp),
    }
    return pd.DataFrame(cols)


def reproduce_fig2(tp: pd.DataFrame, surv: pd.DataFrame) -> pd.DataFrame:
    """KM by MGMT under several cohort/OS definitions. The paper does not say which it used."""
    s = surv.merge(tp[["id", "mgmt"]], left_on="cohort_tp_id", right_on="id", how="left")
    s["os_event_provided"] = 1  # authors' OS: treat every provided value as a death
    # Unrepaired v3: OS sits on other rows than v1's, so take v3's own OS rows.
    v3 = tp[tp["v3_os_days_provided"].notna()].sort_values(["patient_id", "tp"])
    v3 = v3.groupby("patient_id", as_index=False).head(1).assign(v3_os_event=1)
    variants = {
        "provided_subset__v1_os": (
            s[s["in_provided_subset"]],
            "os_days_provided",
            "os_event_provided",
        ),
        "v3_unrepaired__v3_os": (v3, "v3_os_days_provided", "v3_os_event"),
        "rule_cohort__clean_os": (s, "os_days", "os_event"),
        "rule_cohort__recomputed_os": (s, "os_days_recomputed", "os_event"),
    }
    # All patients, OS from dates, MGMT from first timepoint.
    pt_mgmt = tp.sort_values("tp").groupby("patient_id")["mgmt"].first()
    allp = surv.attrs.get("all_patients")
    if allp is not None:
        allp = allp.assign(mgmt=allp["patient_id"].map(pt_mgmt))
        variants["all_patients__recomputed_os"] = (allp, "os_days_recomputed", "os_event")

    rows = []
    for name, (d, time, event) in variants.items():
        d = d[d["mgmt"].isin(["Methylated", "Unmethylated"])]
        p = km_by_group(
            d,
            time,
            event,
            "mgmt",
            f"OS by MGMT — {name}",
            FIGURES_DIR / f"fig2_km_mgmt__{name}.png",
        )
        n = d.dropna(subset=[time])
        rows.append(
            {
                "variant": name,
                "n": len(n),
                "events": int(n[event].sum()),
                "n_meth": int((n["mgmt"] == "Methylated").sum()),
                "n_unmeth": int((n["mgmt"] == "Unmethylated").sum()),
                "logrank_p": p,
                "paper_p": PAPER_FIG2_P,
            }
        )
    return pd.DataFrame(rows)


@app.command()
def main():
    from ucsd_pfs.survival import recompute_os
    from ucsd_pfs.timeline import patient_timeline

    tp = pd.read_parquet(PROCESSED_TIMEPOINTS)
    surv = pd.read_parquet(PROCESSED_SURVIVAL)
    surv.attrs["all_patients"] = recompute_os(patient_timeline(tp))

    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    t1 = reproduce_table1(tp, surv)
    t1.to_csv(TABLES_DIR / "table1_reproduction.csv")
    logger.info("Table 1:\n" + t1.to_string())

    f2 = reproduce_fig2(tp, surv)
    f2.to_csv(TABLES_DIR / "fig2_reproduction.csv", index=False)
    logger.info("Fig. 2 (KM OS by MGMT):\n" + f2.to_string(index=False))


if __name__ == "__main__":
    app()
