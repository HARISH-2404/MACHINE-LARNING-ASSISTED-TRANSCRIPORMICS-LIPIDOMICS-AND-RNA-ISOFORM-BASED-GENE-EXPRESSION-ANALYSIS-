#!/usr/bin/env python3
"""Integrate gene, transcript and lipid matrices; run exploratory ML classification."""
import argparse
from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import make_scorer, balanced_accuracy_score
from sklearn.inspection import permutation_importance

def load_block(path, prefix, sample_ids, log_counts=False):
    df = pd.read_csv(path).set_index("feature_id")
    df = df.loc[:, sample_ids].apply(pd.to_numeric, errors="coerce")
    if log_counts:
        # log transform is for exploratory ML features only; original counts remain unchanged on disk.
        df = np.log2(df + 1)
    df.index = [f"{prefix}::{x}" for x in df.index]
    return df.T

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--metadata", required=True)
    p.add_argument("--gene", required=True)
    p.add_argument("--transcript", required=True)
    p.add_argument("--lipid", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--top-features", type=int, default=100)
    a = p.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(a.metadata).set_index("sample_id")
    if meta.index.duplicated().any():
        raise ValueError("sample_id must be unique.")
    y_raw = meta["condition"].astype(str)
    if y_raw.nunique() != 2:
        raise ValueError("This demo classifier expects exactly two condition classes.")
    samples = list(meta.index)
    blocks = [
        load_block(a.gene, "GENE", samples, True),
        load_block(a.transcript, "TRANSCRIPT", samples, True),
        load_block(a.lipid, "LIPID", samples, True),
    ]
    X = pd.concat(blocks, axis=1)
    # Remove features with no variation; any imputation/scaling/selection occurs inside CV pipeline.
    X = X.loc[:, X.nunique(dropna=True) > 1]
    X = X.replace([np.inf, -np.inf], np.nan)
    if X.shape[1] == 0:
        raise ValueError("No variable features found.")
    enc = LabelEncoder()
    y = enc.fit_transform(y_raw.loc[X.index])
    counts = pd.Series(y).value_counts()
    n_splits = min(5, int(counts.min()))
    if n_splits < 2:
        raise ValueError("At least two samples per class are needed for cross-validation.")
    k = min(a.top_features, X.shape[1])
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=a.seed)
    models = {
        "RandomForest": RandomForestClassifier(n_estimators=400, class_weight="balanced", random_state=a.seed),
        "SVM_RBF": SVC(kernel="rbf", class_weight="balanced", probability=True, random_state=a.seed),
        "ElasticNet_Logistic": LogisticRegression(penalty="elasticnet", solver="saga", l1_ratio=0.5,
                                                  class_weight="balanced", max_iter=8000, random_state=a.seed)
    }
    records = []
    for name, model in models.items():
        pipe = Pipeline([
            ("imputer", SimpleImputer(strategy="median", add_indicator=True)),
            ("scale", StandardScaler()),
            ("select", SelectKBest(score_func=f_classif, k=k)),
            ("model", model)
        ])
        scoring = {"balanced_accuracy": "balanced_accuracy", "accuracy": "accuracy", "roc_auc": "roc_auc"}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            scores = cross_validate(pipe, X, y, cv=cv, scoring=scoring, error_score=np.nan)
        records.append({"model": name, "n_splits": n_splits,
                        "balanced_accuracy_mean": np.nanmean(scores["test_balanced_accuracy"]),
                        "balanced_accuracy_sd": np.nanstd(scores["test_balanced_accuracy"], ddof=1),
                        "accuracy_mean": np.nanmean(scores["test_accuracy"]),
                        "roc_auc_mean": np.nanmean(scores["test_roc_auc"])})
    metrics = pd.DataFrame(records).sort_values("balanced_accuracy_mean", ascending=False)
    metrics.to_csv(out / "model_cv_metrics.csv", index=False)
    # Fit a descriptive final model on all samples for feature ranking only.
    # Its in-sample feature ranks are exploratory and are not an unbiased performance estimate.
    best_name = metrics.iloc[0]["model"]
    pipe = Pipeline([
        ("imputer", SimpleImputer(strategy="median", add_indicator=True)),
        ("scale", StandardScaler()),
        ("select", SelectKBest(score_func=f_classif, k=k)),
        ("model", models[best_name])
    ])
    pipe.fit(X, y)
    selected_mask = pipe.named_steps["select"].get_support()
    selected_names = X.columns[selected_mask]
    model = pipe.named_steps["model"]
    if hasattr(model, "feature_importances_"):
        importance = model.feature_importances_
    elif hasattr(model, "coef_"):
        importance = np.mean(np.abs(model.coef_), axis=0)
    else:
        importance = np.zeros(len(selected_names))
    rank = pd.DataFrame({"feature": selected_names, "importance": importance}).sort_values("importance", ascending=False)
    rank.to_csv(out / "exploratory_feature_importance.csv", index=False)
    plt.figure(figsize=(9, 6))
    plotdf = metrics.sort_values("balanced_accuracy_mean")
    plt.barh(plotdf.model, plotdf.balanced_accuracy_mean, xerr=plotdf.balanced_accuracy_sd.fillna(0))
    plt.xlim(0, 1); plt.xlabel("Mean cross-validated balanced accuracy")
    plt.title("Exploratory model comparison"); plt.tight_layout()
    plt.savefig(out / "model_comparison.png", dpi=200); plt.close()
    plt.figure(figsize=(9, 7))
    top = rank.head(25).sort_values("importance")
    if len(top):
        plt.barh(top.feature, top.importance)
    plt.xlabel("Model-specific feature importance"); plt.title(f"Exploratory top features: {best_name}")
    plt.tight_layout(); plt.savefig(out / "top_features.png", dpi=200); plt.close()
    with open(out / "README_RESULTS.txt", "w", encoding="utf-8") as f:
        f.write("Exploratory results only. Cross-validation estimates are uncertain for small datasets.\n")
        f.write("Feature importance from the final all-sample fit is descriptive, not an unbiased biomarker validation.\n")
        f.write("Best model by mean balanced accuracy: " + best_name + "\n")
        f.write("Classes: " + ", ".join(enc.classes_) + "\n")
    print(f"Integrated {X.shape[0]} samples and {X.shape[1]} variable features.")
    print(metrics.to_string(index=False))
    print(f"Results saved to {out.resolve()}")

if __name__ == "__main__":
    main()
