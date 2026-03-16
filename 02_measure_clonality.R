## ============================================================================
## 02_measure_clonality.R
## --------------------------------------------------------------------------
## Purpose : Correct observed variant allele frequencies (VAF) for tumour
##           content (purity) and local copy number, then derive a per-variant
##           clonality score on a 0-1 scale.
##
## Rationale:
##   Raw VAFs are diluted by admixed normal DNA and distorted by copy-number
##   aberrations.  This script applies a standard correction formula to obtain
##   the "true" cancer-cell frequency of each variant, then doubles it to
##   estimate clonality (under the assumption that most somatic SNVs are
##   heterozygous in a diploid background: a fully clonal het has VAF ~ 0.5).
##
## Input files:
##   - Input/sequenza.rds               (Sequenza CN segments)
##   - Input/SNVS.Rdata                 (somatic SNV calls with observed VAF)
##   - Input/NCG_CRC_cancer_genes.Rdata (CRC driver gene list from NCG)
##
## Output files:
##   - Results/SNVS.clonality.Rdata     (annotated SNV table, R binary)
##   - Results/Mutations.csv            (same table, CSV export)
##   - Results/Frequency.pdf            (raw vs TC-corrected VAF histograms)
##   - Results/Clonality.pdf            (TC-corrected VAF vs clonality histograms)
##   - Results/DensityPlot.pdf          (density of clonality with driver gene labels)
## ============================================================================

library(plyr)      # ddply for split-apply-combine
library(ggplot2)   # density plot and theming
library(ggrepel)   # non-overlapping gene-name labels

# --- Load data ---------------------------------------------------------------
sequenza <- readRDS("Input/sequenza.rds")                # CN segments
load("Input/SNVS.Rdata")                                 # provides: snvs
load("Input/NCG_CRC_cancer_genes.Rdata")                 # provides: crc_cancer_genes

# --- Tumour-content (TC) correction function ---------------------------------
# Formula adjusts the observed VAF for:
#   - tumour content / cellularity (tc): fraction of tumour cells in the sample
#   - local total copy number in the tumour (CNt)
#   - copy number in matched normal tissue (CNn, default = 2 for diploid)
#
# Derivation:
#   observed VAF = (tc * VAF_true * CNt) / (tc * CNt + (1 - tc) * CNn)
#   solving for VAF_true and capping at 1 to avoid artefactually high values.
get.tc.correction.somatic <- function(obs, tc, CNt, CNn = 2) {
  return(min(obs * (1 + ((CNn * (1 - tc)) / (CNt * tc))), 1))
}

# --- Apply correction to every variant ---------------------------------------
for (i in 1:nrow(snvs)) {
  # Compute TC-corrected VAF; total CN = major + minor allele CN at the locus
  snvs$frq.tc[i] <- get.tc.correction.somatic(
    snvs$obs.VAF[i],                             # observed VAF (0-1)
    snvs$Aberrant_Cell_Fraction,                  # sample-level tumour purity
    snvs$major_cn[i] + snvs$minor_cn[i],          # local total CN in tumour
    snvs$normal_cn[i]                              # local CN in normal (usually 2)
  )
  # Clonality = 2 * corrected VAF, capped at 1.
  # A fully clonal heterozygous variant in a diploid tumour has corrected
  # VAF ~ 0.5, so doubling maps it to clonality ~ 1.0.
  snvs$clonality[i] <- min(c(snvs$frq.tc[i] * 2, 1))
}

# --- Persist results ---------------------------------------------------------
save(snvs, file = "Results/SNVS.clonality.Rdata")
write.csv(snvs, file = "Results/Mutations.csv", row.names = FALSE)

# =============================================================================
# PLOTS
# =============================================================================

# --- Histogram pair 1: raw VAF vs TC-corrected VAF ---------------------------
# Shows the effect of purity/CN correction on the frequency distribution.
pdf(file = "Results/Frequency.pdf", h = 5, w = 12)
par(mfrow = c(1, 2))
hist(snvs$obs.VAF, main = "Frequency", xlab = "Frequency",
     ylab = "Counts", xlim = c(0, 1), ylim = c(0, 700))
hist(snvs$frq.tc, main = "Frequency\npost TC correction", xlab = "Frequency",
     ylab = "Counts", xlim = c(0, 1), ylim = c(0, 700))
dev.off()

# --- Histogram pair 2: TC-corrected VAF vs clonality -------------------------
# Illustrates the doubling step that converts corrected VAF to clonality.
pdf(file = "Results/Clonality.pdf", h = 5, w = 12)
par(mfrow = c(1, 2))
hist(snvs$frq.tc, main = "Frequency\npost TC correction", xlab = "Frequency",
     ylab = "Counts", xlim = c(0, 1), ylim = c(0, 700))
hist(snvs$clonality, main = "Clonality", xlab = "Frequency",
     ylab = "Counts", xlim = c(0, 1), ylim = c(0, 700))
dev.off()

# --- Helper: custom x-axis with optional reverse -----------------------------
# Draws a baseline segment and sets breaks/labels; optionally reverses the axis
# so that high clonality (clonal) appears on the left.
base_breaks_x <- function(x, br, la, reverse = TRUE) {
  d <- data.frame(y = -Inf, yend = -Inf, x = 0, xend = 100)
  if (reverse) {
    return(list(geom_segment(data = d, aes(x = x, y = y, xend = xend, yend = yend),
                             inherit.aes = FALSE),
                scale_x_reverse(breaks = br, labels = la)))
  } else {
    return(list(geom_segment(data = d, aes(x = x, y = y, xend = xend, yend = yend),
                             inherit.aes = FALSE),
                scale_x_continuous(breaks = br, labels = la)))
  }
}

# --- Helper: custom y-axis ---------------------------------------------------
# Draws a baseline segment on the right side and sets y-axis breaks/labels.
base_breaks_y <- function(m, br, la) {
  d <- data.frame(x = Inf, xend = Inf, y = 0, yend = m)
  list(geom_segment(data = d, aes(x = x, y = y, xend = xend, yend = yend),
                    inherit.aes = FALSE),
       scale_y_continuous(breaks = br, labels = la))
}

# --- Density plot: clonality distribution with driver-gene annotations -------
# The density is drawn using counts (not probability), so the y-axis reflects
# the expected number of alterations at each clonality level.
# Vertical dashed lines at 0.80 and 0.35 mark the thresholds used in
# 03_clone_composition.R to classify variants as monoclonal / biclonal / polyclonal.
p <- ggplot(snvs, aes(x = clonality)) +
  geom_density(fill = rgb(172, 152, 199, maxColorValue = 255), alpha = 0.5,
               aes(y = after_stat(count))) +
  geom_vline(xintercept = 0.80, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0.35, linetype = "dashed", color = "grey")

# Extract the internally computed density curve (x, y coordinates) so we can
# project driver-gene points onto the curve at matching clonality values.
gg   <- ggplot_build(p)$data[[1]]
ymax <- max(gg$y)
ymed <- ymax / 2
lmax <- ifelse(ymax >= 0.5, round(ymax), 1)
lmed <- ifelse(lmax == 1, 0.5, ifelse(ymed > 0.5, round(ymed), 0.5))

# Apply theme and reverse x-axis (clonal on the left, subclonal on the right)
p <- p + theme_bw() + scale_x_reverse() +
  ylab("Expected Alterations") + xlab("Alteration clonality (%)")

# --- Identify known CRC driver genes among the called variants ---------------
# For each unique clonality value that harbours a driver, collect gene names.
drvrs <- ddply(
  subset(snvs, Hugo_Symbol %in% crc_cancer_genes$symbol),
  .(clonality), summarise, gname = unique(Hugo_Symbol)
)

# Map each driver gene to the nearest (x, y) point on the density curve
# so the label sits directly on the curve.
drvrs$x <- NA
drvrs$y <- NA
for (i in 1:nrow(drvrs)) {
  x   <- drvrs$clonality[i]
  idx <- which(abs(gg$x - x) == min(abs(gg$x - x)))   # closest density point
  drvrs$x[i] <- gg$x[idx]
  drvrs$y[i] <- gg$y[idx]
}

# --- Render the final density plot with gene labels --------------------------
pdf(file = "Results/DensityPlot.pdf", h = 6, w = 9)
p +
  geom_point(data = drvrs, aes(x = x, y = y), colour = "black",
             show.legend = FALSE) +
  geom_text_repel(data = drvrs, aes(x = x, y = y, label = gname))
dev.off()
