#!/usr/bin/env python3
"""Exploratory transcript/isoform usage analysis from transcript abundance matrix.
Input: transcript_counts.csv with feature_id and sample columns; metadata.csv with sample_id, condition.
This computes per-gene isoform proportions and tests condition differences on proportions (Welch t-test).
It is exploratory; use IsoformSwitchAnalyzeR/SUPPA2 for validated switching and consequence analysis.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--transcripts',required=True); p.add_argument('--metadata',required=True)
    p.add_argument('--mapping',required=True,help='CSV with transcript_id,gene_id columns')
    p.add_argument('--outdir',required=True); p.add_argument('--pseudocount',type=float,default=1e-8)
    a=p.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    md=pd.read_csv(a.metadata).set_index('sample_id')
    if 'condition' not in md: raise ValueError('metadata requires condition column')
    if md.index.duplicated().any(): raise ValueError('sample_id must be unique')
    groups=list(md.condition.astype(str).unique())
    if len(groups)!=2: raise ValueError('Exactly two conditions are required for this exploratory test')
    df=pd.read_csv(a.transcripts).set_index('feature_id')
    mapping=pd.read_csv(a.mapping)
    if not {'transcript_id','gene_id'}.issubset(mapping.columns): raise ValueError('mapping requires transcript_id,gene_id')
    samples=list(md.index)
    missing=[s for s in samples if s not in df.columns]
    if missing: raise ValueError(f'Metadata samples missing from transcript matrix: {missing[:5]}')
    mat=df[samples].apply(pd.to_numeric,errors='coerce').fillna(0).clip(lower=0)
    mat.index=mat.index.astype(str)
    mapping=mapping.drop_duplicates('transcript_id').set_index('transcript_id')
    mapping=mapping.loc[mapping.index.intersection(mat.index)]
    mat=mat.loc[mapping.index]
    props=pd.DataFrame(index=mat.index,columns=samples,dtype=float)
    for gene, ids in mapping.groupby('gene_id').groups.items():
        block=mat.loc[list(ids)]
        denom=block.sum(axis=0)+a.pseudocount
        props.loc[list(ids)]=block.div(denom,axis=1)
    props.index.name='transcript_id'; props.to_csv(out/'isoform_usage_proportions.csv')
    c1=md.index[md.condition.astype(str)==groups[0]]; c2=md.index[md.condition.astype(str)==groups[1]]
    rows=[]
    for tid,row in props.iterrows():
        x=row[c1].astype(float).values; y=row[c2].astype(float).values
        stat,pv=ttest_ind(x,y,equal_var=False,nan_policy='omit')
        rows.append({'transcript_id':tid,'gene_id':mapping.loc[tid,'gene_id'],'mean_usage_'+groups[0]:np.nanmean(x),'mean_usage_'+groups[1]:np.nanmean(y),'delta_usage':np.nanmean(y)-np.nanmean(x),'t_statistic':stat,'p_value':pv})
    res=pd.DataFrame(rows); res['padj_BH']=multipletests(res.p_value.fillna(1),method='fdr_bh')[1]
    res.sort_values('padj_BH').to_csv(out/'isoform_usage_tests_exploratory.csv',index=False)
    top=res.nsmallest(min(25,len(res)),'padj_BH').sort_values('delta_usage')
    plt.figure(figsize=(9,max(4,len(top)*.25))); plt.barh(top.transcript_id,top.delta_usage); plt.axvline(0,linewidth=.8)
    plt.xlabel('Change in mean isoform proportion ('+groups[1]+' − '+groups[0]+')'); plt.title('Exploratory transcript usage differences')
    plt.tight_layout(); plt.savefig(out/'top_isoform_usage_changes.png',dpi=200); plt.close()
    (out/'INTERPRETATION.txt').write_text('Exploratory proportion-level Welch tests; not a substitute for replicate-aware differential transcript usage / isoform switching. Confirm transcript-to-gene mapping and quantification method. Use IsoformSwitchAnalyzeR or SUPPA2 for formal switching and functional consequences.\n',encoding='utf-8')
    print(f'Wrote exploratory usage results for {len(res)} transcripts to {out}')
if __name__=='__main__': main()
