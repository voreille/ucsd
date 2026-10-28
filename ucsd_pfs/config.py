"""Project-wide paths and constants."""

from pathlib import Path

from loguru import logger

PROJ_ROOT = Path(__file__).resolve().parents[1]

DATA_DIR = PROJ_ROOT / "data"
RAW_DATA_DIR = DATA_DIR / "raw"
INTERIM_DATA_DIR = DATA_DIR / "interim"
PROCESSED_DATA_DIR = DATA_DIR / "processed"
EXTERNAL_DATA_DIR = DATA_DIR / "external"

MODELS_DIR = PROJ_ROOT / "models"
REPORTS_DIR = PROJ_ROOT / "reports"
FIGURES_DIR = REPORTS_DIR / "figures"
TABLES_DIR = REPORTS_DIR / "tables"

TCIA_UPLOADS = "https://www.cancerimagingarchive.net/wp-content/uploads"
CLINICAL_FILES = {
    "v1": "UCSD_PTGBM-clinical-information_v1_2025-25-03.xlsx",
    "v3": "UCSD_PTGBM-clinical-information_v3_2026-12-Mar.xlsx",
}
NEGATIVE_CASES_FILE = "details_of_negative_cases_TCIA.xlsx"
CLINICAL_SHEET = "Cinical info"  # sic, typo is in the source file

INTERIM_TIMEPOINTS = INTERIM_DATA_DIR / "clinical_timepoints.parquet"
PROCESSED_TIMEPOINTS = PROCESSED_DATA_DIR / "timepoints.parquet"
PROCESSED_SURVIVAL = PROCESSED_DATA_DIR / "survival_patients.parquet"

# Paper: "OS and PFS information is available for a subset of 94 patients with an MRI scan
# with RSI within 90 days following surgery."
POSTOP_WINDOW_DAYS = 90
# Tolerance (days) when comparing provided OS with OS recomputed from date offsets.
OS_MATCH_TOL_DAYS = 3

SEED = 42

logger.debug(f"PROJ_ROOT is {PROJ_ROOT}")
