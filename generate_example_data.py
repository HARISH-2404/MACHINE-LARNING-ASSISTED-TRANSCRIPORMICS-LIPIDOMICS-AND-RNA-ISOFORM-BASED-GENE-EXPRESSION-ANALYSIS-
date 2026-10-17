#!/usr/bin/env python3
"""Generate deterministic synthetic multi-omics data for software testing only."""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "example"
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(42)

n_per_group = 12
samples = [f"C{i:02d}" for i in range(1, n_per_group+1)] + [f"S{i:02d}" for i in range(1, n_per_group+1)]
conditions = ["Control"] * n_per_group + ["Stress"] * n_per_group
batches = [f"B{1 + (i % 2)}" for i in range(len(samples))]
meta = pd.DataFrame({"sample_id": samples, "condition": conditions, "batch": batches})
meta.to_csv(OUT / "metadata.csv", index=False)

def counts_table(prefix, n_features, effect_n):
    base = rng.uniform(30, 450, size=n_features)
    mat = np.zeros((n_features, len(samples)), dtype=int)
    for j in range(len(samples)):
        mu = base.copy()
        if j >= n_per_group:
            mu[:effect_n] *= 3.0
            mu[effect_n:2*effect_n] *= 0.35
        # Negative-binomial-like overdispersed counts via gamma-Poisson mixture
        lam = rng.gamma(shape=8.0, scale=mu/8.0)
        mat[:, j] = rng.poisson(lam)
    ids = [f"{prefix}{i:03d}" for i in range(1, n_features+1)]
    return pd.DataFrame(mat, columns=samples).assign(feature_id=ids).loc[:, ["feature_id"] + samples]

counts_table("Gene", 250, 35).to_csv(OUT / "gene_counts.csv", index=False)
counts_table("TX", 300, 30).to_csv(OUT / "transcript_counts.csv", index=False)

lipids = []
for i in range(1, 81):
    base = rng.lognormal(mean=10, sigma=0.5, size=len(samples))
    if i <= 12:
        base[n_per_group:] *= 2.0
    elif i <= 24:
        base[n_per_group:] *= 0.5
    noise = rng.lognormal(mean=0, sigma=0.15, size=len(samples))
    lipids.append(base * noise)
lipid_df = pd.DataFrame(lipids, columns=samples)
lipid_df.insert(0, "feature_id", [f"Lipid_{i:03d}" for i in range(1, 81)])
lipid_df.to_csv(OUT / "lipidomics.csv", index=False)
print(f"Synthetic test data written to: {OUT}")
