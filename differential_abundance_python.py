#!/usr/bin/env python3
"""Exploratory Welch tests for gene/transcript/lipid feature matrices with BH-FDR.
For raw RNA-seq integer counts use DESeq2/edgeR instead; this script is not a replacement.
"""
import argparse
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

def analyse(matrix_path, metadata, outdir, block_name, log_transform=False):
    df=pd.read_csv(matrix_path).set_index('feature_id')
    samples=list(metadata.index)
    miss=[x for x in samples if x not in df.columns]
    if miss: raise ValueError(f'{matrix_path}: missing sample columns {miss[:5]}')
    x=df[samples].apply(pd.to_numeric,errors='coerce')
    if log_transform: x=np.log2(x.clip(lower=0)+1)
    groups=list(metadata.condition.astype(str).unique())
    if len(groups)!=2: raise ValueError('Exactly two conditions required')
    g1=metadata.index[metadata.condition.astype(str)==groups[0]]; g2=metadata.index[metadata.condition.astype(str)==groups[1]]
    rows=[]
    for fid,row in x.iterrows():
        a=row[g1].dropna().astype(float).values; b=row[g2].dropna().astype(float).values
        if len(a)<2 or len(b)<2: pv=np.nan; stat=np.nan
        else: stat,pv=ttest_ind(a,b,equal_var=False,nan_policy='omit')
        rows.append({'feature_id':fid,'mean_'+groups[0]:np.mean(a) if len(a) else np.nan,'mean_'+groups[1]:np.mean(b) if len(b) else np.nan,'difference_'+groups[1]+'_minus_'+groups[0]:(np.mean(b)-np.mean(a)) if len(a) and len(b) else np.nan,'t_statistic':stat,'p_value':pv})
    out=pd.DataFrame(rows); valid=out.p_value.notna(); out['padj_BH']=np.nan
    if valid.any(): out.loc[valid,'padj_BH']=multipletests(out.loc[valid,'p_value'],method='fdr_bh')[1]
    out=out.sort_values('padj_BH',na_position='last'); out.to_csv(outdir/f'{block_name}_exploratory_differential_abundance.csv',index=False)

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--metadata',required=True); p.add_argument('--gene'); p.add_argument('--transcript'); p.add_argument('--lipid'); p.add_argument('--outdir',required=True); a=p.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True); md=pd.read_csv(a.metadata).set_index('sample_id')
    if 'condition' not in md: raise ValueError('metadata requires condition')
    for name,path,log in [('gene',a.gene,True),('transcript',a.transcript,True),('lipid',a.lipid,True)]:
        if path: analyse(path,md,out,name,log_transform=log)
    (out/'README.txt').write_text('These are exploratory Welch tests on log2(feature+1) values with Benjamini-Hochberg correction. For raw RNA-seq counts, use DESeq2 or edgeR. Ensure biological replicates and appropriate design.\n',encoding='utf-8')
    print('Finished exploratory differential abundance.')
if __name__=='__main__': main()
