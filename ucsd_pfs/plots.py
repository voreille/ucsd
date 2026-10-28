"""Plotting helpers."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from lifelines import KaplanMeierFitter  # noqa: E402
from lifelines.statistics import logrank_test  # noqa: E402


def km_by_group(
    df: pd.DataFrame,
    time: str,
    event: str,
    group: str,
    title: str,
    out: Path | None = None,
    time_unit_days: float = 30.44,
) -> float:
    """Kaplan-Meier curves for a two-level grouping; returns the log-rank p-value."""
    d = df.dropna(subset=[time, event, group])
    levels = sorted(d[group].unique())
    if len(levels) != 2:
        raise ValueError(f"Expected 2 groups in {group}, got {levels}")
    a, b = (d[d[group] == lv] for lv in levels)
    p = logrank_test(a[time], b[time], a[event], b[event]).p_value

    fig, ax = plt.subplots(figsize=(6, 4.5))
    for lv, sub in zip(levels, (a, b), strict=True):
        km = KaplanMeierFitter(label=f"{lv} (n={len(sub)})")
        km.fit(sub[time] / time_unit_days, sub[event])
        km.plot_survival_function(ax=ax, ci_show=True)
    ax.set_xlabel("Months since initial surgery")
    ax.set_ylabel("Survival probability")
    ax.set_title(f"{title}\nlog-rank p = {p:.2g}")
    ax.set_ylim(0, 1.02)
    fig.tight_layout()
    if out is not None:
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=150)
    plt.close(fig)
    return float(p)
