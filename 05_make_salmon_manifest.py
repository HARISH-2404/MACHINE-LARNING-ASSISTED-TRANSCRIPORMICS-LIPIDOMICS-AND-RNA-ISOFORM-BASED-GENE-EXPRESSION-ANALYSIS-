#!/usr/bin/env python3
"""Validate/standardize a sample manifest for the Salmon shell workflow."""
import argparse,csv
from pathlib import Path
p=argparse.ArgumentParser(); p.add_argument('--input',required=True); p.add_argument('--output',default='data/input_templates/salmon_manifest.csv'); a=p.parse_args()
rows=list(csv.DictReader(open(a.input,newline='',encoding='utf-8-sig')))
need={'sample_id','fastq_r1','fastq_r2'}
if not rows or not need.issubset(rows[0]): raise SystemExit('Need columns: sample_id,fastq_r1,fastq_r2')
ids=[r['sample_id'].strip() for r in rows]
if any(not x for x in ids) or len(set(ids))!=len(ids): raise SystemExit('sample_id values must be nonempty and unique')
for r in rows:
 for k in ['fastq_r1','fastq_r2']:
  f=Path(r[k]); f=f if f.is_absolute() else Path(a.input).resolve().parent/f
  if not f.exists(): raise SystemExit(f'Missing {k} for {r["sample_id"]}: {f}')
out=Path(a.output); out.parent.mkdir(parents=True,exist_ok=True)
with out.open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=['sample_id','fastq_r1','fastq_r2']); w.writeheader(); w.writerows(rows)
print('Validated manifest written to',out)
