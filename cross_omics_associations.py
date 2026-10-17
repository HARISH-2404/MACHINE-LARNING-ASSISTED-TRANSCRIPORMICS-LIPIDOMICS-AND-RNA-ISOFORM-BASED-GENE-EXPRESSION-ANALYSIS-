#!/usr/bin/env python3
"""Exploratory cross-omics feature associations using matched sample IDs and FDR correction."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

def load(path, prefix):
    df=pd.read_csv(path).set_index("feature_id").apply(pd.to_numeric, errors="coerce")
    df.index=[f"{prefix}::{i}" for i in df.index]
    return df

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--gene',required=True); p.add_argument('--transcript',required=True); p.add_argument('--lipid',required=True)
    p.add_argument('--metadata',required=True); p.add_argument('--outdir',required=True)
    p.add_argument('--max-features',type=int,default=50); p.add_argument('--min-abs-rho',type=float,default=0.5)
    a=p.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    meta=pd.read_csv(a.metadata).set_index('sample_id')
    blocks=[load(a.gene,'GENE'),load(a.transcript,'TX'),load(a.lipid,'LIPID')]
    samples=list(meta.index)
    for b in blocks:
        if not set(samples).issubset(b.columns): raise ValueError('Sample IDs do not match across metadata and omics data')
    # Use transformed feature values for exploratory association, not inferential DE testing.
    gene=np.log2(blocks[0][samples]+1); tx=np.log2(blocks[1][samples]+1); lipid=np.log2(blocks[2][samples]+1)
    # Limit to highest-variance features to bound pairwise computation; document this selection.
    def topvar(df): return df.loc[df.var(axis=1).sort_values(ascending=False).head(a.max_features).index]
    gene,tx,lipid=map(topvar,[gene,tx,lipid]); results=[]
    for block_name,left,right in [('gene_lipid',gene,lipid),('transcript_lipid',tx,lipid),('gene_transcript',gene,tx)]:
        for fid,x in left.iterrows():
            for gid,y in right.iterrows():
                ok=x.notna() & y.notna()
                if ok.sum()<5 or x[ok].nunique()<2 or y[ok].nunique()<2: continue
                rho,pv=spearmanr(x[ok],y[ok])
                results.append((block_name,fid,gid,rho,pv,int(ok.sum())))
    res=pd.DataFrame(results,columns=['comparison','feature_a','feature_b','spearman_rho','pvalue','n_samples'])
    if len(res):
        res['padj']=multipletests(res.pvalue,method='fdr_bh')[1]
        res=res.sort_values('padj'); res.to_csv(out/'cross_omics_correlations.csv',index=False)
        selected=res[(res.padj<0.05)&(res.spearman_rho.abs()>=a.min_abs_rho)]
        selected.to_csv(out/'candidate_cross_omics_associations.csv',index=False)
        top=res.head(30).sort_values('spearman_rho')
        plt.figure(figsize=(10,7)); plt.barh(top.feature_a+' ↔ '+top.feature_b,top.spearman_rho)
        plt.axvline(0,linewidth=0.8); plt.xlabel('Spearman correlation'); plt.title('Top cross-omics associations (exploratory)')
        plt.tight_layout(); plt.savefig(out/'cross_omics_top_correlations.png',dpi=200); plt.close()
    else:
        pd.DataFrame(columns=['comparison','feature_a','feature_b','spearman_rho','pvalue','n_samples','padj']).to_csv(out/'cross_omics_correlations.csv',index=False)
    print(f'Wrote exploratory cross-omics associations to {out.resolve()}')
if __name__=='__main__': main()
