#!/usr/bin/env python3
# RF STRICT MODELS — ACCESSORY GENES (1–99% FILTERED)
# Chromosome | Plasmid | Combined
# Prints CV metrics + Saves everything
# ============================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import joblib
from pathlib import Path

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.metrics import (
    accuracy_score, roc_auc_score, precision_score,
    recall_score, f1_score, confusion_matrix, roc_curve
)

# ============================================================
# BASE PATH
# ============================================================

BASE_DIR = Path(__file__).resolve().parents[2]
PROCESSED_DIR = BASE_DIR / "processed"
FS_DIR = PROCESSED_DIR / "feature_selection/results"

# ============================================================
# LOAD FILTERED MATRICES
# ============================================================

chrom_df_filtered = pd.read_csv(PROCESSED_DIR / "chromosome_ML_filtered_1-99.csv")
plasmid_df_filtered = pd.read_csv(PROCESSED_DIR / "plasmid_ML_filtered_1-99.csv")
combined_df_filtered = pd.read_csv(PROCESSED_DIR / "combined_ML_filtered_1-99.csv")


# ============================================================
# RF TRAINING FUNCTION
# ============================================================

def run_rf_strict(df_filtered, dataset_name):

    print("\n===================================================")
    print(f"Running RF STRICT — {dataset_name.upper()}")
    print("===================================================")

    OUT_DIR = BASE_DIR / f"results/Models/RF/{dataset_name}/STRICT"
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------------
    # Load STRICT feature list
    # --------------------------------------------------------
    strict_features = pd.read_csv(
        FS_DIR / f"{dataset_name}_FINAL_STRICT_features.csv"
    )["gene"].tolist()

    pd.DataFrame({"gene": strict_features}).to_csv(
        OUT_DIR / f"{dataset_name}_STRICT_feature_list_used.csv",
        index=False
    )

    # --------------------------------------------------------
    # Define X and y
    # --------------------------------------------------------
    X = df_filtered[strict_features]
    y = df_filtered["phenotype"].astype(int)

    # --------------------------------------------------------
    # Train/Test Split
    # --------------------------------------------------------
    X_train, X_test, y_train, y_test = train_test_split(
        X, y,
        test_size=0.2,
        stratify=y,
        random_state=42
    )

    # --------------------------------------------------------
    # Model
    # --------------------------------------------------------
    rf = RandomForestClassifier(
        n_estimators=500,
        n_jobs=-1,
        random_state=42,
        class_weight="balanced"
    )

    rf.fit(X_train, y_train)

    joblib.dump(rf, OUT_DIR / f"RF_{dataset_name}_STRICT_model.pkl")

    # --------------------------------------------------------
    # TEST METRICS
    # --------------------------------------------------------
    y_pred = rf.predict(X_test)
    y_prob = rf.predict_proba(X_test)[:, 1]

    test_acc  = accuracy_score(y_test, y_pred)
    test_auc  = roc_auc_score(y_test, y_prob)
    test_prec = precision_score(y_test, y_pred)
    test_rec  = recall_score(y_test, y_pred)
    test_f1   = f1_score(y_test, y_pred)
    cm        = confusion_matrix(y_test, y_pred)

    print("\n--- TEST SET PERFORMANCE ---")
    print(f"Accuracy : {test_acc:.4f}")
    print(f"AUC      : {test_auc:.4f}")
    print(f"Precision: {test_prec:.4f}")
    print(f"Recall   : {test_rec:.4f}")
    print(f"F1-score : {test_f1:.4f}")

    # --------------------------------------------------------
    # 5-FOLD CV ON TRAIN SET
    # --------------------------------------------------------
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    cv_acc, cv_auc, cv_prec, cv_rec, cv_f1 = [], [], [], [], []

    for tr, va in cv.split(X_train, y_train):

        rf.fit(X_train.iloc[tr], y_train.iloc[tr])

        yp  = rf.predict(X_train.iloc[va])
        ypr = rf.predict_proba(X_train.iloc[va])[:, 1]

        cv_acc.append(accuracy_score(y_train.iloc[va], yp))
        cv_auc.append(roc_auc_score(y_train.iloc[va], ypr))
        cv_prec.append(precision_score(y_train.iloc[va], yp))
        cv_rec.append(recall_score(y_train.iloc[va], yp))
        cv_f1.append(f1_score(y_train.iloc[va], yp))

    print("\n--- 5-FOLD CV PERFORMANCE (TRAIN SET) ---")
    print(f"Accuracy : {np.mean(cv_acc):.4f} ± {np.std(cv_acc):.4f}")
    print(f"AUC      : {np.mean(cv_auc):.4f} ± {np.std(cv_auc):.4f}")
    print(f"Precision: {np.mean(cv_prec):.4f} ± {np.std(cv_prec):.4f}")
    print(f"Recall   : {np.mean(cv_rec):.4f} ± {np.std(cv_rec):.4f}")
    print(f"F1-score : {np.mean(cv_f1):.4f} ± {np.std(cv_f1):.4f}")

    # --------------------------------------------------------
    # SAVE METRICS
    # --------------------------------------------------------
    pd.DataFrame([{
        "Model": "RF",
        "Dataset": dataset_name,
        "Feature_Set": "STRICT",
        "Test_Accuracy": test_acc,
        "Test_AUC": test_auc,
        "Test_Precision": test_prec,
        "Test_Recall": test_rec,
        "Test_F1": test_f1,
        "CV_Accuracy_Mean": np.mean(cv_acc),
        "CV_Accuracy_STD": np.std(cv_acc),
        "CV_AUC_Mean": np.mean(cv_auc),
        "CV_AUC_STD": np.std(cv_auc),
        "CV_Precision_Mean": np.mean(cv_prec),
        "CV_Precision_STD": np.std(cv_prec),
        "CV_Recall_Mean": np.mean(cv_rec),
        "CV_Recall_STD": np.std(cv_rec),
        "CV_F1_Mean": np.mean(cv_f1),
        "CV_F1_STD": np.std(cv_f1)
    }]).to_csv(
        OUT_DIR / f"RF_{dataset_name}_STRICT_metrics.csv",
        index=False
    )

    # --------------------------------------------------------
    # CONFUSION MATRIX
    # --------------------------------------------------------
    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues")
    plt.title(f"{dataset_name} RF STRICT — Confusion Matrix")
    plt.savefig(OUT_DIR / f"RF_{dataset_name}_STRICT_CM.png", dpi=300)
    plt.close()

    # --------------------------------------------------------
    # ROC CURVE
    # --------------------------------------------------------
    fpr, tpr, _ = roc_curve(y_test, y_prob)
    plt.plot(fpr, tpr, label=f"AUC={test_auc:.3f}")
    plt.plot([0, 1], [0, 1], 'k--')
    plt.title(f"{dataset_name} RF STRICT — ROC")
    plt.legend()
    plt.savefig(OUT_DIR / f"RF_{dataset_name}_STRICT_ROC.png", dpi=300)
    plt.close()

    print(f"\n✔ {dataset_name.upper()} STRICT model completed successfully")


# ============================================================
# RUN ALL THREE DATASETS
# ============================================================

run_rf_strict(chrom_df_filtered, "chromosome")
run_rf_strict(plasmid_df_filtered, "plasmid")
run_rf_strict(combined_df_filtered, "combined")

print("\n🎉 All RF STRICT models trained and saved successfully.")
# ============================================================
# Random Forest — RELAXED MODELS
# Chromosome | Plasmid | Combined
# Accessory genes (1–99% filtered)
# ============================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import joblib
from pathlib import Path

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.metrics import (
    accuracy_score, roc_auc_score, precision_score,
    recall_score, f1_score, confusion_matrix, roc_curve
)

# ============================================================
# PATHS
# ============================================================

BASE_DIR = Path(__file__).resolve().parents[2]
PROCESSED_DIR = BASE_DIR / "processed"
FS_DIR = PROCESSED_DIR / "feature_selection/results"

# ============================================================
# FUNCTION
# ============================================================

def run_rf_relaxed(dataset_name):

    print("\n===================================================")
    print(f"Running RF RELAXED — {dataset_name.upper()}")
    print("===================================================")

    OUT_DIR = BASE_DIR / f"results/Models/RF/{dataset_name}/RELAXED"
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------------
    # Load Data
    # --------------------------------------------------------
    df = pd.read_csv(
        PROCESSED_DIR / f"{dataset_name}_ML_filtered_1-99.csv"
    )

    relaxed_features = pd.read_csv(
        FS_DIR / f"{dataset_name}_FINAL_RELAXED_features.csv"
    )["gene"].tolist()

    # Save feature list used
    pd.DataFrame({"gene": relaxed_features}).to_csv(
        OUT_DIR / f"RF_{dataset_name}_RELAXED_feature_list_used.csv",
        index=False
    )

    X = df[relaxed_features]
    y = df["phenotype"].astype(int)

    # --------------------------------------------------------
    # Train/Test Split
    # --------------------------------------------------------
    X_train, X_test, y_train, y_test = train_test_split(
        X, y,
        test_size=0.2,
        stratify=y,
        random_state=42
    )

    # --------------------------------------------------------
    # Model
    # --------------------------------------------------------
    rf = RandomForestClassifier(
        n_estimators=500,
        n_jobs=-1,
        random_state=42,
        class_weight="balanced"
    )

    rf.fit(X_train, y_train)

    # Save model
    joblib.dump(
        rf,
        OUT_DIR / f"RF_{dataset_name}_RELAXED_model.pkl"
    )

    # ========================================================
    # TEST METRICS
    # ========================================================

    y_pred = rf.predict(X_test)
    y_prob = rf.predict_proba(X_test)[:, 1]

    test_acc  = accuracy_score(y_test, y_pred)
    test_auc  = roc_auc_score(y_test, y_prob)
    test_prec = precision_score(y_test, y_pred)
    test_rec  = recall_score(y_test, y_pred)
    test_f1   = f1_score(y_test, y_pred)

    print("\n--- TEST SET PERFORMANCE ---")
    print(f"Accuracy : {test_acc:.4f}")
    print(f"AUC      : {test_auc:.4f}")
    print(f"Precision: {test_prec:.4f}")
    print(f"Recall   : {test_rec:.4f}")
    print(f"F1-score : {test_f1:.4f}")

    # ========================================================
    # 5-FOLD CV
    # ========================================================

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    cv_acc, cv_auc, cv_prec, cv_rec, cv_f1 = [], [], [], [], []

    for tr, va in cv.split(X_train, y_train):

        rf.fit(X_train.iloc[tr], y_train.iloc[tr])

        yp  = rf.predict(X_train.iloc[va])
        ypr = rf.predict_proba(X_train.iloc[va])[:, 1]

        cv_acc.append(accuracy_score(y_train.iloc[va], yp))
        cv_auc.append(roc_auc_score(y_train.iloc[va], ypr))
        cv_prec.append(precision_score(y_train.iloc[va], yp))
        cv_rec.append(recall_score(y_train.iloc[va], yp))
        cv_f1.append(f1_score(y_train.iloc[va], yp))

    print("\n--- 5-FOLD CV PERFORMANCE (TRAIN SET) ---")
    print(f"Accuracy : {np.mean(cv_acc):.4f} ± {np.std(cv_acc):.4f}")
    print(f"AUC      : {np.mean(cv_auc):.4f} ± {np.std(cv_auc):.4f}")
    print(f"Precision: {np.mean(cv_prec):.4f} ± {np.std(cv_prec):.4f}")
    print(f"Recall   : {np.mean(cv_rec):.4f} ± {np.std(cv_rec):.4f}")
    print(f"F1-score : {np.mean(cv_f1):.4f} ± {np.std(cv_f1):.4f}")

    # ========================================================
    # SAVE METRICS CSV
    # ========================================================

    metrics_df = pd.DataFrame([{
        "Model": "RF",
        "Dataset": dataset_name,
        "Feature_Set": "RELAXED",
        "Test_Accuracy": test_acc,
        "Test_AUC": test_auc,
        "Test_Precision": test_prec,
        "Test_Recall": test_rec,
        "Test_F1": test_f1,
        "CV_Accuracy_Mean": np.mean(cv_acc),
        "CV_Accuracy_STD": np.std(cv_acc),
        "CV_AUC_Mean": np.mean(cv_auc),
        "CV_AUC_STD": np.std(cv_auc),
        "CV_Precision_Mean": np.mean(cv_prec),
        "CV_Precision_STD": np.std(cv_prec),
        "CV_Recall_Mean": np.mean(cv_rec),
        "CV_Recall_STD": np.std(cv_rec),
        "CV_F1_Mean": np.mean(cv_f1),
        "CV_F1_STD": np.std(cv_f1)
    }])

    metrics_df.to_csv(
        OUT_DIR / f"RF_{dataset_name}_RELAXED_metrics.csv",
        index=False
    )

    # ========================================================
    # SAVE METRICS IMAGE (SAFE — NO ROUND ERROR)
    # ========================================================

    metrics_display = metrics_df.copy()
    numeric_cols = metrics_display.select_dtypes(include=[np.number]).columns
    metrics_display[numeric_cols] = metrics_display[numeric_cols].round(4)

    plt.figure(figsize=(14, 4))
    plt.axis('off')

    table = plt.table(
        cellText=metrics_display.astype(str).values,
        colLabels=metrics_display.columns,
        loc='center'
    )

    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.5)

    plt.savefig(
        OUT_DIR / f"RF_{dataset_name}_RELAXED_metrics.png",
        dpi=300,
        bbox_inches='tight'
    )
    plt.close()

    # ========================================================
    # CONFUSION MATRIX
    # ========================================================

    cm = confusion_matrix(y_test, y_pred)

    plt.figure()
    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues")
    plt.title(f"{dataset_name} RF RELAXED — Confusion Matrix")
    plt.savefig(
        OUT_DIR / f"RF_{dataset_name}_RELAXED_CM.png",
        dpi=300
    )
    plt.close()

    # ========================================================
    # ROC CURVE
    # ========================================================

    fpr, tpr, _ = roc_curve(y_test, y_prob)

    plt.figure()
    plt.plot(fpr, tpr, label=f"AUC={test_auc:.3f}")
    plt.plot([0, 1], [0, 1], 'k--')
    plt.title(f"{dataset_name} RF RELAXED — ROC")
    plt.legend()

    plt.savefig(
        OUT_DIR / f"RF_{dataset_name}_RELAXED_ROC.png",
        dpi=300
    )
    plt.close()

    print(f"\n✔ {dataset_name.upper()} RELAXED model completed successfully")


# ============================================================
# RUN ALL THREE
# ============================================================

for ds in ["chromosome", "plasmid", "combined"]:
    run_rf_relaxed(ds)

print("\n🎉 All RF RELAXED models trained and saved successfully.")
