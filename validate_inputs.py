#!/usr/bin/env python3
"""Validate sample IDs and numeric matrix formats before analysis."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

def read_matrix(path, name, integer=False):
    df = pd.read_csv(path)
    if df.shape[1] < 3:
        raise ValueError(f"{name}: expected feature_id plus at least two sample columns.")
    if df.columns[0] != "feature_id":
        raise ValueError(f"{name}: first column must be named feature_id.")
    if df["feature_id"].isna().any() or df["feature_id"].duplicated().any():
        raise ValueError(f"{name}: feature_id values must be present and unique.")
    values = df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
    if values.isna().all(axis=None):
        raise ValueError(f"{name}: no numeric values found.")
    if (values < 0).any(axis=None):
        raise ValueError(f"{name}: negative values found.")
    if integer and (values.isna().any(axis=None) or not np.allclose(values.values, np.round(values.values))):
        raise ValueError(f"{name}: gene/transcript counts must be non-missing integers.")
    return df, values

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--data-dir", default="data/example")
    a = p.parse_args()
    d = Path(a.data_dir)
    meta = pd.read_csv(d / "metadata.csv")
    required = {"sample_id", "condition"}
    if not required.issubset(meta.columns):
        raise ValueError("metadata.csv must contain sample_id and condition.")
    if meta.sample_id.isna().any() or meta.sample_id.duplicated().any():
        raise ValueError("metadata sample_id values must be present and unique.")
    if meta.condition.nunique() != 2:
        raise ValueError("This demo expects exactly two condition groups.")
    sample_ids = set(meta.sample_id.astype(str))
    for fname, integer in [("gene_counts.csv", True), ("transcript_counts.csv", True), ("lipidomics.csv", False)]:
        df, vals = read_matrix(d / fname, fname, integer)
        cols = set(df.columns[1:].astype(str))
        if cols != sample_ids:
            raise ValueError(f"{fname}: sample columns do not match metadata IDs. Missing={sample_ids-cols}; extra={cols-sample_ids}")
        print(f"OK {fname}: {df.shape[0]} features x {vals.shape[1]} samples")
    print(f"OK metadata.csv: {len(meta)} samples; groups={meta.condition.value_counts().to_dict()}")
    print("Validation passed.")

if __name__ == "__main__":
    main()
