#!/usr/bin/env python3
"""Create a concise Markdown summary from generated analysis result tables."""
import argparse
from pathlib import Path
import pandas as pd

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--results',required=True); p.add_argument('--out',default='analysis_summary.md'); a=p.parse_args()
    root=Path(a.results); lines=['# Multi-omics analysis summary','', '> This report summarizes files found in the supplied results directory. Verify methods, design, and biological interpretation before publication.','']
    for path in sorted(root.rglob('*.csv')):
        try: df=pd.read_csv(path)
        except Exception: continue
        lines += [f'## {path.relative_to(root)}','',f'- Rows: {len(df)}','- Columns: '+', '.join(map(str,df.columns))]
        pcol=next((c for c in ['padj','FDR','qvalue'] if c in df.columns),None)
        if pcol:
            n=int((pd.to_numeric(df[pcol],errors='coerce')<0.05).sum()); lines.append(f'- Features with {pcol} < 0.05: {n}')
        lines.append('')
    out=Path(a.out); out.parent.mkdir(parents=True,exist_ok=True); out.write_text('\n'.join(lines),encoding='utf-8')
    print(f'Report written: {out.resolve()}')
if __name__=='__main__': main()
