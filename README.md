# Lung Adenocarcinoma RNA-seq: From Raw Reads to Candidate Biomarkers

An end-to-end, reproducible RNA-seq analysis of lung adenocarcinoma (LUAD), built in two parts:

1. **Sequencing pipeline (Linux / Snakemake):** raw paired-end FASTQ → QC (FastQC, fastp, MultiQC) → transcript quantification (Salmon) → gene-level counts (tximport). Runs locally or on an HPC cluster via a SLURM profile.
2. **Cancer genomics analysis (R / Python / SQL):** TCGA-LUAD tumor vs. matched normal differential expression (DESeq2, paired design), pathway enrichment (fgsea, MSigDB Hallmarks), survival association of top genes (Cox / Kaplan–Meier), results stored in a SQLite database and served through an interactive R Shiny explorer.

> **Question:** Which genes and pathways most consistently distinguish LUAD tumors from adjacent normal lung tissue in the same patient, and are any of them associated with overall survival?

---

## Repository layout

```
├── config/                 # Snakemake config + sample sheet
├── workflow/Snakefile      # FASTQ → QC → Salmon pipeline
├── profiles/slurm/         # Run the pipeline on an HPC cluster
├── scripts/
│   ├── 00_fetch_fastq.sh   # Download + subsample public FASTQ (SRA)
│   ├── 01_tximport.R       # Salmon transcripts → gene counts
│   ├── 10_download_tcga.R  # TCGA-LUAD STAR counts via TCGAbiolinks
│   ├── 11_deseq2_paired.R  # Paired tumor vs. normal DE + figures
│   ├── 12_enrichment.R     # Hallmark GSEA
│   ├── 13_survival.R       # Cox / KM for top DE genes
│   └── 20_build_db.py      # Load results into SQLite
├── app/app.R               # Shiny explorer backed by SQLite
├── environment.yml         # Conda environment (pipeline tools)
└── results/                # Figures and tables (generated)
```

## How to run

```bash
# 1. Environment
conda env create -f environment.yml
conda activate luad-rnaseq

# 2. Sequencing pipeline (fill in config/samples.tsv first)
bash scripts/00_fetch_fastq.sh
snakemake --cores 4                      # local
snakemake --profile profiles/slurm       # HPC cluster
Rscript scripts/01_tximport.R

# 3. TCGA analysis
Rscript scripts/10_download_tcga.R
Rscript scripts/11_deseq2_paired.R
Rscript scripts/12_enrichment.R
Rscript scripts/13_survival.R
python scripts/20_build_db.py

# 4. Explorer
Rscript -e 'shiny::runApp("app")'
```

## Key results

_Fill in after running — 3–5 bullets with real numbers, each pointing to a figure._

- Paired DESeq2 on **N** matched tumor/normal pairs identified **X** genes at FDR < 0.05 and |log2FC| > 1 (`results/figures/volcano.png`).
- Samples separate cleanly by tissue type on PC1 (**Y%** of variance) (`results/figures/pca.png`).
- Top enriched Hallmark pathways: … (`results/figures/gsea_hallmarks.png`).
- Of the top 25 DE genes, **Z** were associated with overall survival after FDR correction (`results/figures/km_<gene>.png`).

## Methods notes and limitations

- Paired design (`~ patient + condition`) controls for between-patient variation; only patients with both a primary tumor and a solid-tissue-normal sample are included.
- Survival analysis uses a median expression split for visualization and continuous Cox models for testing; results are exploratory and not validated in an independent cohort.
- The FASTQ pipeline is run on a subsampled public dataset to keep compute small; the same workflow scales to full data on a cluster.

## Skills demonstrated

Linux/Bash · Snakemake · SLURM · Conda · Git · FastQC/fastp/MultiQC · Salmon · R (DESeq2, tximport, fgsea, survival, ggplot2) · Python (pandas, sqlite3) · SQL · R Shiny
