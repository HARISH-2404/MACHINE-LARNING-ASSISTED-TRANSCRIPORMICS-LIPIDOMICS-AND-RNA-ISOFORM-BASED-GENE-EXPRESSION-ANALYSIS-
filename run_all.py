#!/usr/bin/env python3
"""Run the Python demonstration workflow from repository root."""
import subprocess, sys
from pathlib import Path
root = Path(__file__).resolve().parents[1]
def run(args):
    print("\n>>>", " ".join(map(str, args)), flush=True)
    subprocess.run(list(map(str, args)), cwd=root, check=True)
py = sys.executable
run([py, "scripts/generate_example_data.py"])
run([py, "scripts/validate_inputs.py", "--data-dir", "data/example"])
run([py, "scripts/lipidomics_analysis.py", "--data", "data/example/lipidomics.csv",
     "--metadata", "data/example/metadata.csv", "--outdir", "results/lipidomics"])
run([py, "scripts/integrate_ml.py", "--metadata", "data/example/metadata.csv",
     "--gene", "data/example/gene_counts.csv", "--transcript", "data/example/transcript_counts.csv",
     "--lipid", "data/example/lipidomics.csv", "--outdir", "results/integration"])
print("\nPython demonstration completed. Run the R scripts separately for DESeq2 analyses.")
