# Clonal Evolution Tutorial

A step-by-step R pipeline for analysing **clonal evolution** in a colorectal cancer (CRC) whole-exome sequencing sample.  The tutorial walks through copy-number profiling, variant-allele-frequency correction, clonality estimation, clone-composition classification, and VAF-based mutation clustering.

The example data comes from **TCGA sample TCGA-CA-6718** (primary CRC, matched blood normal).

---

## Pipeline overview

| Script | Step | Description |
|--------|------|-------------|
| `00_configuration.R` | Setup | Install all required R packages (CRAN, Bioconductor, GitHub) |
| `01_analyse_CN_profile.R` | CN profiling | Visualise allele-specific copy-number segments from Sequenza |
| `02_measure_clonality.R` | Clonality scoring | Correct VAF for tumour content and CN; derive clonality (0–1) |
| `03_clone_composition.R` | Clone classification | Classify variants as monoclonal / biclonal / polyclonal |
| `04_clustering_mutations.R` | Mutation clustering | Cluster variants into subclones using sciClone |

Scripts are designed to be executed **sequentially** (00 → 01 → 02 → 03 → 04), as each step depends on outputs from the previous one.

---

## Repository structure

```
ClonalEvolution_tutorial/
├── 00_configuration.R
├── 01_analyse_CN_profile.R
├── 02_measure_clonality.R
├── 03_clone_composition.R
├── 04_clustering_mutations.R
├── Input/                        # Pre-computed input data
│   ├── sequenza.rds              # Sequenza CN segment calls
│   ├── sequenza.Rdata            # Sequenza results (R binary)
│   ├── ASCAT.Rdata               # ASCAT CN analysis results
│   ├── SNVS.Rdata                # Somatic SNV calls with VAF and CN annotation
│   ├── NCG_CRC_cancer_genes.Rdata# CRC driver gene list (Network of Cancer Genes)
│   └── sciClone.Rdata            # Pre-computed sciClone clustering (fallback)
├── Results/                      # Pipeline outputs (generated at runtime)
└── Help_output/                  # Reference outputs for validation
    ├── sequenza.png              # Expected CN profile plot
    ├── sequenza.results.csv      # Expected CN segment table
    ├── ASCAT.png                 # ASCAT CN profile (reference)
    ├── ASCAT.results.csv         # ASCAT segment table (reference)
    ├── CNV.genes.csv             # CN-altered cancer genes
    ├── Mutations.csv             # Annotated mutation table with clonality
    ├── SNVS.clonality.Rdata      # Clonality-annotated SNVs (R binary)
    ├── Frequency.pdf             # Raw vs TC-corrected VAF histograms
    ├── Clonality.pdf             # Corrected VAF vs clonality histograms
    ├── DensityPlot.pdf           # Clonality density with driver gene labels
    ├── CloneComposition.pdf      # Stacked bar chart of clone fractions
    ├── clones.xls                # sciClone cluster assignments
    └── Clones.stats.csv          # Per-cluster summary + driver events
```

---

## Dependencies

| Package | Source | Purpose |
|---------|--------|---------|
| `plyr` | CRAN | Split-apply-combine operations (`ddply`) |
| `ggplot2` | CRAN | Plotting (density, bar charts) |
| `reshape2` | CRAN | Wide ↔ long data reshaping (`melt`) |
| `ggrepel` | CRAN | Non-overlapping text labels on plots |
| `RColorBrewer` | CRAN | Colour palettes |
| `devtools` | CRAN | GitHub package installation |
| `IRanges` | Bioconductor | Interval/range operations |
| `bmm` | GitHub (`genome/bmm`) | Bayesian mixture models (sciClone dependency) |
| `sciClone` | GitHub (`genome/sciClone`) | VAF-based subclone clustering |

Run `00_configuration.R` to install all dependencies automatically.

---

## Step-by-step description

### 00 — Configuration

Installs CRAN packages, the Bioconductor `IRanges` package (via `BiocManager`), and GitHub packages `bmm` and `sciClone` (via `devtools`).

### 01 — Copy-number profile

Reads the Sequenza segment table and produces:

- **`sequenza.results.csv`** — full segment-level CN table (chromosome, start, end, major/minor allele CN, ploidy, cellularity).
- **`sequenza.png`** — genome-wide plot with major-allele CN in blue and minor-allele CN in yellow.  Segments exceeding CN = 5 are drawn in a dimmed colour.

### 02 — Clonality measurement

Corrects each variant's observed VAF for tumour content (purity) and local copy number using the formula:

```
VAF_corrected = VAF_obs × [ 1 + CNn × (1 − tc) / (CNt × tc) ]
```

where `tc` = tumour content, `CNt` = local tumour CN, `CNn` = normal CN (default 2).  Clonality is then defined as `min(2 × VAF_corrected, 1)`, under the assumption that most somatic SNVs are heterozygous in a diploid background.

Outputs include:
- Annotated mutation table (`Mutations.csv` / `SNVS.clonality.Rdata`)
- Histogram comparisons (raw vs corrected VAF; corrected VAF vs clonality)
- Density plot of clonality with CRC driver genes (e.g. *KRAS*, *PIK3CA*, *EGFR*) projected onto the curve

### 03 — Clone composition

Classifies each variant into three categories based on clonality thresholds:

| Category | Clonality range | Interpretation |
|----------|----------------|----------------|
| Monoclonal | ≥ 0.80 | Present in virtually all tumour cells |
| Biclonal | 0.35 – 0.80 | Present in a major subpopulation |
| Polyclonal | < 0.35 | Present in a minor subpopulation |

Produces a horizontal stacked bar chart (`CloneComposition.pdf`) showing the proportion of variants in each category per patient.

### 04 — Mutation clustering

Uses **sciClone** to cluster variants by VAF into discrete subclones:

- Variants in CN-aberrant regions are excluded (their VAF is distorted by gains/losses).
- A Bayesian mixture model is fitted to the remaining VAFs (up to 4 clusters, minimum depth = 50 reads).
- Each cluster is annotated with summary statistics (median clonality, range) and known driver events for this sample.

---

## Usage

```r
# 1. Install dependencies (run once)
source("00_configuration.R")

# 2. Run the analysis pipeline
source("01_analyse_CN_profile.R")
source("02_measure_clonality.R")
source("03_clone_composition.R")
source("04_clustering_mutations.R")
```

Set your working directory to the repository root before running:

```r
setwd("/path/to/ClonalEvolution_tutorial")
```

All outputs are written to the `Results/` directory.  Compare them against the reference files in `Help_output/` to verify correctness.

---

## References

- **Sequenza** — Favero F et al. *Sequenza: allele-specific copy number and mutation profiles from tumor sequencing data.* Ann Oncol. 2015;26(1):64-70.
- **sciClone** — Miller CA et al. *SciClone: inferring clonal architecture and tracking the spatial and temporal patterns of tumor evolution.* PLoS Comput Biol. 2014;10(8):e1003665.
- **NCG** — Network of Cancer Genes: https://network-cancer-genes.org/
