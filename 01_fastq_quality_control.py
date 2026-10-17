#!/usr/bin/env python3
"""Run FastQC, optional fastp trimming, and MultiQC on paired-end FASTQ files.
Input manifest CSV columns: sample_id, fastq_r1, fastq_r2. Paths may be relative to manifest.
"""
import argparse, csv, shutil, subprocess, sys
from pathlib import Path

def run(cmd):
    print("+", " ".join(map(str, cmd)), flush=True)
    subprocess.run(list(map(str, cmd)), check=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--outdir', default='results/fastq_qc')
    ap.add_argument('--trim', action='store_true', help='Run fastp before final FastQC')
    ap.add_argument('--threads', type=int, default=4)
    a=ap.parse_args(); out=Path(a.outdir).resolve(); out.mkdir(parents=True,exist_ok=True)
    for tool in ['fastqc','multiqc']:
        if not shutil.which(tool): raise SystemExit(f'Required executable not found on PATH: {tool}')
    rows=list(csv.DictReader(open(a.manifest,newline='',encoding='utf-8-sig')))
    if not rows or not {'sample_id','fastq_r1','fastq_r2'} <= set(rows[0]):
        raise SystemExit('Manifest must have sample_id,fastq_r1,fastq_r2 columns.')
    qc_inputs=[]; trimmed=out/'trimmed'; trimmed.mkdir(exist_ok=True)
    for r in rows:
        s=r['sample_id'].strip()
        r1=Path(r['fastq_r1']); r2=Path(r['fastq_r2'])
        if not r1.is_absolute(): r1=(Path(a.manifest).resolve().parent/r1).resolve()
        if not r2.is_absolute(): r2=(Path(a.manifest).resolve().parent/r2).resolve()
        if not r1.exists() or not r2.exists(): raise SystemExit(f'Missing FASTQ for {s}: {r1} / {r2}')
        if a.trim:
            if not shutil.which('fastp'): raise SystemExit('--trim requested but fastp not found on PATH')
            o1=trimmed/f'{s}_R1.trimmed.fastq.gz'; o2=trimmed/f'{s}_R2.trimmed.fastq.gz'
            run(['fastp','-i',r1,'-I',r2,'-o',o1,'-O',o2,'--thread',a.threads,'--html',out/f'{s}.fastp.html','--json',out/f'{s}.fastp.json'])
            r1,r2=o1,o2
        qc_inputs += [r1,r2]
    run(['fastqc','--threads',a.threads,'--outdir',out,*qc_inputs])
    run(['multiqc',out,'--outdir',out/'multiqc','--force'])
    print(f'FASTQ QC completed. Review {out}/multiqc/multiqc_report.html')
if __name__=='__main__': main()
