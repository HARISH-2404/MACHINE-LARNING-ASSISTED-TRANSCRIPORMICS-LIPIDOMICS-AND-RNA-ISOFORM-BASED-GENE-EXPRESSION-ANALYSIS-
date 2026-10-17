#!/usr/bin/env python3
"""Basic lipidomics QC, median sample normalization, log2 transform, Welch tests and FDR."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--reference", default=None, help="Reference condition; defaults to first metadata level.")
    p.add_argument("--comparison", default=None, help="Comparison condition; defaults to second metadata level.")
    a = p.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    x = pd.read_csv(a.data).set_index("feature_id")
    meta = pd.read_csv(a.metadata).set_index("sample_id")
    if not set(meta.index).issubset(x.columns):
        raise ValueError("Metadata sample IDs are missing from lipidomics matrix.")
    x = x.loc[:, meta.index].apply(pd.to_numeric, errors="coerce")
    groups = list(meta["condition"].dropna().unique())
    if len(groups) != 2:
        raise ValueError("Exactly two condition groups are required.")
    ref = a.reference or groups[0]; comp = a.comparison or groups[1]
    if ref not in groups or comp not in groups or ref == comp:
        raise ValueError(f"Invalid comparison. Available groups: {groups}")
    # Remove features missing in >50% of all samples. No silent zero substitution.
    x = x.loc[x.isna().mean(axis=1) <= 0.5]
    # For this baseline demo, impute remaining missing values with half the smallest positive value per feature.
    # For real studies, choose a method appropriate to the missingness mechanism.
    for idx in x.index:
        row = x.loc[idx]
        positive = row[row > 0]
        fill = positive.min() / 2 if len(positive) else np.nan
        x.loc[idx] = row.fillna(fill)
    x = x.dropna(axis=0, how="any")
    # Median sample normalization. Inspect QC and choose normalization appropriate to study.
    sample_medians = x.median(axis=0)
    if (sample_medians <= 0).any():
        raise ValueError("Non-positive sample median after missing-value handling.")
    norm = x.div(sample_medians, axis=1).mul(float(sample_medians.median()))
    logx = np.log2(norm + 1)
    logx.to_csv(out / "lipidomics_log2_normalized.csv")
    rows = []
    for feature, row in logx.iterrows():
        v1 = row[meta.index[meta.condition == ref]].dropna()
        v2 = row[meta.index[meta.condition == comp]].dropna()
        if len(v1) < 2 or len(v2) < 2:
            continue
        stat, pval = ttest_ind(v2, v1, equal_var=False)
        rows.append({"feature_id": feature, "mean_reference": v1.mean(),
                     "mean_comparison": v2.mean(), "log2FC_comparison_vs_reference": v2.mean()-v1.mean(),
                     "pvalue": pval})
    res = pd.DataFrame(rows)
    if len(res):
        res["padj"] = multipletests(res.pvalue, method="fdr_bh")[1]
        res = res.sort_values("padj")
    res.to_csv(out / "lipidomics_differential_results.csv", index=False)
    # Heatmap of up to 30 most variable features
    top = logx.var(axis=1).sort_values(ascending=False).head(min(30, len(logx))).index
    if len(top):
        plt.figure(figsize=(max(9, len(meta)*0.4), max(5, len(top)*0.25)))
        sns.heatmap(logx.loc[top], cmap="vlag", xticklabels=True, yticklabels=True)
        plt.title("Lipidomics: most variable features (log2 normalized)")
        plt.tight_layout(); plt.savefig(out / "lipidomics_heatmap.png", dpi=200); plt.close()
    plt.figure(figsize=(7, 5))
    for group in [ref, comp]:
        vals = logx.loc[:, meta.index[meta.condition == group]].mean(axis=1)
        plt.hist(vals, bins=20, alpha=0.5, label=group)
    plt.xlabel("Mean log2 normalized abundance per lipid"); plt.ylabel("Number of lipid features")
    plt.legend(); plt.tight_layout(); plt.savefig(out / "lipidomics_distribution.png", dpi=200); plt.close()
    print(f"Results saved to {out.resolve()}")

if __name__ == "__main__":
    main()
