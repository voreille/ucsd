#################################################################################
# GLOBALS                                                                       #
#################################################################################

PROJECT_NAME = ucsd_pfs
PYTHON = uv run python

#################################################################################
# COMMANDS                                                                      #
#################################################################################

.PHONY: requirements data clean-survival reproduce notebooks all test lint format clean

## Install dependencies
requirements:
	uv sync

## Download clinical spreadsheets and build data/interim/clinical_timepoints.parquet
data:
	$(PYTHON) -m $(PROJECT_NAME).dataset

## Build cleaned OS/PFS targets in data/processed/
clean-survival: data
	$(PYTHON) -m $(PROJECT_NAME).survival

## Reproduce paper Table 1 and Fig. 2 into reports/
reproduce: clean-survival
	$(PYTHON) -m $(PROJECT_NAME).reproduce

## Execute all notebooks in place
notebooks: clean-survival
	uv run jupyter nbconvert --to notebook --execute --inplace notebooks/*.ipynb

## Full pipeline
all: reproduce notebooks

## Run unit tests
test:
	uv run pytest

## Lint
lint:
	uv run ruff check .
	uv run ruff format --check .

## Format
format:
	uv run ruff check --fix .
	uv run ruff format .

## Delete compiled Python files and generated outputs (keeps raw data)
clean:
	find . -type d -name "__pycache__" -prune -exec rm -rf {} +
	rm -f data/interim/*.parquet data/processed/*.parquet reports/tables/*.csv reports/figures/*.png
