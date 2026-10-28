"""Train PFS models on clinical features from the first post-surgery timepoint.

Not implemented yet: survival targets must first be settled (see references/survival_cleaning.md).
Planned: lifelines CoxPH/KM, scikit-survival (CoxNet, RSF, GBSA) and sklearn baselines, evaluated
with repeated stratified K-fold (N is ~90) and reported with confidence intervals.
"""

import typer

app = typer.Typer()


@app.command()
def main():
    raise NotImplementedError("PFS modeling is the next step; see CLAUDE.md")


if __name__ == "__main__":
    app()
