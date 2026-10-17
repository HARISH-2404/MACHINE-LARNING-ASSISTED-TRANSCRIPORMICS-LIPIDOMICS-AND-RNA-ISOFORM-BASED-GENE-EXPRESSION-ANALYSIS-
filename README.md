# AI-Driven Multi-Omics Integration Pipeline

**Transcriptomics + RNA isoform profiling + lipidomics + machine learning**

A reproducible starter repository for plant molecular biology research. It includes:
- RNA-seq differential gene expression with DESeq2 (R)
- Transcript/isoform differential abundance using Salmon quantifications and tximport (R)
- Lipidomics quality control, log transformation, differential abundance and plots (Python)
- Cross-omics feature integration and machine-learning classification (Python)
- Synthetic example data generator for software testing

> **Research integrity:** the included data are synthetic and are only for testing. Replace them with real, biologically replicated data before drawing scientific conclusions. The starter pipeline does not automatically prove causal relationships or validate biomarkers.

## Repository structure

```text
AI-MultiOmics-Integration/
├── README.md
├── LICENSE
├── CITATION.cff
├── requirements.txt
├── environment.yml
├── .gitignore
├── config/
│   └── config.yaml
├── data/
│   ├── example/
│   │   ├── metadata.csv
│   │   ├── gene_counts.csv
│   │   ├── transcript_counts.csv
│   │   └── lipidomics.csv
│   └── input_templates/
│       ├── metadata_template.csv
│       ├── gene_counts_template.csv
│       ├── transcript_counts_template.csv
│       └── lipidomics_template.csv
├── scripts/
│   ├── generate_example_data.py
│   ├── validate_inputs.py
│   ├── lipidomics_analysis.py
│   ├── integrate_ml.py
│   └── run_all.py
├── R/
│   ├── differential_gene_expression.R
│   └── transcript_isoform_analysis.R
├── docs/
│   ├── INPUT_FORMATS.md
│   ├── WORKFLOW.md
│   └── LIMITATIONS.md
└── results/
    └── .gitkeep
```

## 1. Requirements

Recommended: Python 3.10–3.12 and R 4.3 or later. The Python scripts may work on other compatible versions, but package compatibility is not guaranteed.

### Python setup

```bash
python -m venv .venv
```

Windows:
```powershell
.venv\Scripts\activate
```

macOS/Linux:
```bash
source .venv/bin/activate
```

Then:
```bash
python -m pip install --upgrade pip
pip install -r requirements.txt
```

### R setup

Open R or RStudio and install Bioconductor packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("DESeq2", "tximport", "tximeta"))
install.packages(c("readr", "dplyr", "tibble"))
```

For isoform switching / alternative splicing, additional tools may be needed, such as `IsoformSwitchAnalyzeR`, `SUPPA2`, `StringTie`, or long-read tools. See `docs/WORKFLOW.md`.

## 2. Run the demonstration

From the repository root:

```bash
python scripts/generate_example_data.py
python scripts/validate_inputs.py --data-dir data/example
python scripts/lipidomics_analysis.py --data data/example/lipidomics.csv --metadata data/example/metadata.csv --outdir results/lipidomics
python scripts/integrate_ml.py --metadata data/example/metadata.csv --gene data/example/gene_counts.csv --transcript data/example/transcript_counts.csv --lipid data/example/lipidomics.csv --outdir results/integration
```

The example has 24 synthetic samples split across two groups, with 12 biological samples per group. This is a code demonstration only; real experimental design and sample size should be justified biologically and statistically.

Run gene-level DESeq2 from R:

```bash
Rscript R/differential_gene_expression.R data/example/gene_counts.csv data/example/metadata.csv results/transcriptomics
```

Run transcript-level exploratory analysis from R:

```bash
Rscript R/transcript_isoform_analysis.R data/example/transcript_counts.csv data/example/metadata.csv results/isoforms
```

The transcript-level script is a simplified transcript-feature DESeq2 analysis. For rigorous transcript quantification, use Salmon/Kallisto with a transcriptome reference and tximport, or a validated long-read pipeline. Transcript differential abundance is not by itself proof of isoform switching or alternative splicing.

## 3. Input design

All datasets must use the same sample IDs as `metadata.csv`. Metadata requires:
- `sample_id`: unique sample identifier
- `condition`: two groups for the included pairwise demo (e.g., Control, Stress)
- `batch`: optional batch variable; use `NA` if unavailable

See `docs/INPUT_FORMATS.md` for complete requirements.

## 4. Scientific interpretation

- DESeq2 uses raw integer read counts; do not supply TPM/FPKM to the gene-count script.
- The lipidomics script uses feature intensities and performs median-based sample normalization followed by log2 transformation. For real studies, choose normalization and missing-value handling based on the instrument, QC samples, and data-generating process.
- The ML script uses repeated cross-validation and reports balanced accuracy, ROC-AUC when possible, and feature importance. A tiny dataset can produce unstable metrics. Keep biological replicates independent and avoid leakage.
- For a real predictive study, do feature selection inside each training fold, use an external test set where possible, and do not use technical replicates as independent biological samples.
- Cross-omics associations are exploratory, not proof of direct molecular regulation.

## 5. License and citation

This project is provided under the MIT License. See `CITATION.cff` for citation metadata. Cite the software versions and reference genome/transcriptome used in any publication.


## Advanced modules included

See `docs/DETAILED_MODULES.md` for detailed instructions. Additional scripts include:

- `scripts/build_salmon_index.sh` — build a Salmon transcriptome index.
- `scripts/run_rnaseq_salmon.sh` — FastQC, optional fastp trimming, Salmon quantification and MultiQC.
- `R/tximport_deseq2.R` — transcript quantification import and gene-level DESeq2.
- `R/isoform_switch_analysis.R` — input checks and preparation for a transcript-switch workflow.
- `scripts/cross_omics_associations.py` — FDR-adjusted cross-omics Spearman associations.
- `scripts/summary_report.py` — builds a Markdown index of result tables.

### Run advanced Python modules

```bash
python scripts/cross_omics_associations.py --gene data/example/gene_counts.csv --transcript data/example/transcript_counts.csv --lipid data/example/lipidomics.csv --metadata data/example/metadata.csv --outdir results/cross_omics
python scripts/summary_report.py --results results --out results/analysis_summary.md
```

For raw paired-end RNA-seq, prepare a tab-separated `samples.tsv` (`sample_id`, `read1`, `read2`), build a Salmon index with `scripts/build_salmon_index.sh`, then run `scripts/run_rnaseq_salmon.sh`. The raw-read workflow requires external command-line tools and a matching transcriptome FASTA.

## 3. Additional Python analyses

### Exploratory differential abundance (gene, transcript, lipid matrices)

```bash
python scripts/differential_abundance_python.py --metadata data/example/metadata.csv --gene data/example/gene_counts.csv --transcript data/example/transcript_counts.csv --lipid data/example/lipidomics.csv --outdir results/python_differential
```

This applies Welch tests to log2(feature + 1) values and Benjamini-Hochberg correction. **Do not use it as the primary differential-expression method for raw RNA-seq counts**; use the DESeq2 R script instead.

### Exploratory isoform usage proportions

```bash
python scripts/isoform_usage_analysis.py --transcripts data/example/transcript_counts.csv --metadata data/example/metadata.csv --mapping data/example/transcript_gene_map.csv --outdir results/isoform_usage
```

The included transcript-to-gene mapping is synthetic for software testing. Replace it with a mapping derived from the exact transcript annotation used for quantification. Proportion-level Welch tests are exploratory and are not equivalent to formal isoform switching.

## 4. Java utilities

Java source code is in `java/`. It provides dependency-free CSV input validation, per-feature summary statistics, and transcript usage proportions. See [`java/README.md`](java/README.md) for compilation and commands. Requires JDK 17 or later.

Compile:

```bash
javac -d java/build java/src/main/java/org/multiomics/MultiOmicsToolkit.java
```

Validate all example inputs:

```bash
java -cp java/build org.multiomics.MultiOmicsToolkit validate data/example/metadata.csv data/example/gene_counts.csv data/example/transcript_counts.csv data/example/lipidomics.csv
```

Feature summary:

```bash
java -cp java/build org.multiomics.MultiOmicsToolkit summarize data/example/lipidomics.csv results/java_lipid_summary.csv
```

Transcript usage proportions:

```bash
java -cp java/build org.multiomics.MultiOmicsToolkit isoform-usage data/example/transcript_counts.csv data/example/transcript_gene_map.csv results/java_isoform_usage.csv
```
