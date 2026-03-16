## ============================================================================
## 04_clustering_mutations.R
## --------------------------------------------------------------------------
## Purpose : Cluster somatic mutations into subclones by their variant allele
##           frequency (VAF) distribution using sciClone, then annotate each
##           cluster with summary statistics and known driver events.
##
## sciClone (Miller et al., PLoS Comput Biol 2014) fits a Bayesian mixture
## model to the VAF distribution while excluding copy-number-aberrant regions,
## yielding discrete mutation clusters that correspond to cell populations
## (clones) present in the tumour at different cellular prevalences.
##
## Input files:
##   - Input/sequenza.rds   (Sequenza CN segments – used to mask CN-altered loci)
##   - Input/SNVS.Rdata     (somatic SNV calls with read counts and observed VAF)
##
## Output files:
##   - Results/clones.xls        (per-variant cluster assignments, tab-delimited)
##   - Results/Clones.stats.csv  (summary statistics per cluster + driver events)
## ============================================================================

library(sciClone)   # VAF-based subclone clustering
library(plyr)       # ddply for per-cluster summaries

# --- Load data ---------------------------------------------------------------
sequenza <- readRDS("Input/sequenza.rds")   # CN segments from Sequenza
load("Input/SNVS.Rdata")                    # provides: snvs

# --- Prepare sciClone input: VAF table ---------------------------------------
# sciClone expects five columns in this order:
#   1. chromosome
#   2. position
#   3. reference-supporting read counts
#   4. variant-supporting read counts
#   5. variant allele fraction on a 0-100 scale
vafs         <- snvs[, c('Chromosome', 'Position', 'ref_counts', 'var_counts', 'obs.VAF')]
vafs$obs.VAF <- vafs$obs.VAF * 100   # convert from 0-1 fraction to 0-100 percentage

# --- Prepare sciClone input: copy-number segments ----------------------------
# sciClone expects four columns:
#   1. chromosome
#   2. segment start position
#   3. segment stop position
#   4. total copy number
# Variants falling in regions with aberrant CN are excluded from clustering
# because their VAF is distorted.  Unlisted regions are assumed diploid (CN = 2).
cnv <- sequenza[, c('Chr', 'Start', 'End', 'cn')]

# --- Run sciClone clustering ------------------------------------------------
# If sciClone is installed, perform de novo clustering; otherwise, load
# pre-computed results shipped with the tutorial.
lib <- installed.packages()
if ("sciClone" %in% rownames(lib)) {
  sc <- sciClone(vafs             = vafs,
                 copyNumberCalls  = cnv,
                 sampleNames      = snvs$Patient,   # sample identifier
                 minimumDepth     = 50,              # ignore variants with <50 reads
                 useSexChrs       = FALSE,           # exclude sex chromosomes
                 maximumClusters  = 4)               # upper bound on cluster count
  # save(sc, file = "Input/sciClone.Rdata")
} else {
  load("Input/sciClone.Rdata")
}

# --- Export cluster assignments and map back to the SNV table ----------------
# writeClusterTable produces a tab-delimited file with cluster IDs per variant.
writeClusterTable(sc, "Results/clones.xls")
tmp <- read.delim("Results/clones.xls")

# Match sciClone output back to the original SNV table using genomic coordinates
# (chromosome + position) as a composite key.
x <- paste0(snvs$Chromosome, ".", snvs$Position)
y <- paste0(tmp[, 1], ".", tmp[, 2])
snvs$cluster <- tmp$cluster[match(x, y)]

# --- Per-cluster summary statistics ------------------------------------------
# For each cluster, compute a six-number summary of clonality values
# (Min, 1st Qu., Median, Mean, 3rd Qu., Max).
stats <- ddply(snvs, .(cluster), function(x) summary(x$clonality))

# --- Annotate clusters with known driver events ------------------------------
# These are sample-specific annotations derived from prior biological knowledge
# of the CRC case under study (TCGA-CA-6718).
stats$events    <- NA
stats$events[1] <- "EGFR amplification"
stats$events[2] <- "PIK3CA chr3:178952086:T>G c=0.72"
stats$events[3] <- "PIK3CA chr3:178916876:G>A c=0.56; KRAS chr12:25398255 G>T c=0.47"

write.csv(stats, "Results/Clones.stats.csv", row.names = FALSE)
