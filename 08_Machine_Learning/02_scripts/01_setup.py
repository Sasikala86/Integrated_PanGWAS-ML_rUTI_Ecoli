#!/usr/bin/env python3
# ============================================================
# ML_Accessory — PROJECT SETUP
# ============================================================
# Creates the repository directories required by the pipeline.
# Raw input files are intentionally not included in this repository.
# ============================================================
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parents[1]

DIRECTORIES = [
    BASE_DIR / "data",
    BASE_DIR / "processed",
    BASE_DIR / "results",
    BASE_DIR / "results" / "EDA",
    BASE_DIR / "results" / "plots",
    BASE_DIR / "processed" / "feature_selection" / "results",
    BASE_DIR / "processed" / "feature_selection" / "plots",
    BASE_DIR / "results" / "Models",
    BASE_DIR / "results" / "SHAP",
]

for directory in DIRECTORIES:
    directory.mkdir(parents=True, exist_ok=True)

print("ML_Accessory directory structure is ready.")
print(f"Repository: {BASE_DIR}")
