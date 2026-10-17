#!/usr/bin/env python3
"""QC plots and summaries for gene/transcript count matrices; does not perform DE testing."""
import argparse
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

def main():
    p=argparse.ArgumentParser(); p.add_argument('--counts',required=True); p.add_argument('--metadata',required=True); p.add_argument('--outdir',default='results/count_qc'); a=p.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    df=pd.read_csv(a.counts); meta=pd.read_csv(a.metadata)
    if 'feature_id' not in df or not {'sample_id','condition'} <= set(meta): raise SystemExit('Counts need feature_id + sample columns; metadata needs sample_id,condition.')
    if meta.sample_id.duplicated().any() or df.feature_id.duplicated().any(): raise SystemExit('Duplicate sample_id or feature_id detected.')
    ids=meta.sample_id.astype(str).tolist(); missing=set(ids)-set(df.columns); extra=set(df.columns[1:])-set(ids)
    if missing or extra: raise SystemExit(f'Sample mismatch. Missing={missing}; extra={extra}')
    x=df.set_index('feature_id')[ids].apply(pd.to_numeric,errors='raise')
    if (x<0).any().any() or not np.allclose(x.values,np.round(x.values)): raise SystemExit('Raw count matrix must contain non-negative integer counts.')
    lib=x.sum(axis=0); summary=pd.DataFrame({'sample_id':ids,'library_size':lib.values,'detected_features':(x>0).sum(axis=0).values}).merge(meta,on='sample_id')
    summary.to_csv(out/'sample_library_qc.csv',index=False)
    plt.figure(figsize=(10,5)); plt.bar(summary.sample_id,summary.library_size); plt.xticks(rotation=90); plt.ylabel('Total raw counts'); plt.tight_layout(); plt.savefig(out/'library_sizes.png',dpi=180); plt.close()
    # Exploratory log-count PCA. This is a QC view, not differential-expression inference.
    mat=np.log2(x.T+1); mat=mat.loc[ids]
    ncomp=min(2,mat.shape[0]-1,mat.shape[1])
    if ncomp>=2:
        scores=PCA(n_components=2).fit_transform(StandardScaler().fit_transform(mat))
        pc=pd.DataFrame(scores,columns=['PC1','PC2']); pc['sample_id']=ids; pc=pc.merge(meta,on='sample_id'); pc.to_csv(out/'pca_scores.csv',index=False)
        plt.figure(figsize=(7,6))
        for group,g in pc.groupby('condition'):
            plt.scatter(g.PC1,g.PC2,label=str(group));
            for _,r in g.iterrows(): plt.annotate(r.sample_id,(r.PC1,r.PC2),fontsize=7)
        plt.xlabel('PC1'); plt.ylabel('PC2'); plt.legend(); plt.tight_layout(); plt.savefig(out/'pca_counts.png',dpi=180); plt.close()
    print(f'Count-matrix QC complete: {out.resolve()}')
if __name__=='__main__': main()
