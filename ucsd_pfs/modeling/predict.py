"""Score a processed feature table with a trained model (not implemented yet)."""

import typer

app = typer.Typer()


@app.command()
def main():
    raise NotImplementedError("PFS modeling is the next step; see CLAUDE.md")


if __name__ == "__main__":
    app()
