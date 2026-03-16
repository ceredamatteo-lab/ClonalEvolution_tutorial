## ============================================================================
## 01_analyse_CN_profile.R
## --------------------------------------------------------------------------
## Purpose : Visualise the allele-specific copy-number (CN) profile produced
##           by Sequenza and export the segment table.
##
## Input files:
##   - Input/sequenza.rds              (Sequenza segment-level CN calls)
##   - Input/NCG_CRC_cancer_genes.Rdata (Network of Cancer Genes – CRC subset)
##
## Output files:
##   - Results/sequenza.results.csv    (segment-level CN table)
##   - Results/sequenza.png            (genome-wide allele-specific CN plot)
##
## The plot displays major-allele CN in blue and minor-allele CN in yellow,
## with chromosome boundaries drawn as vertical grey lines.  Segments whose
## CN exceeds the y-axis limit (5) are drawn in a dimmed colour.
## ============================================================================

library(plyr)   # ddply for per-chromosome probe counts

# --- Load data ---------------------------------------------------------------
# sequenza : data.frame; one row per CN segment.
#   Key columns: Chr, Start, End, nProbes (segment size in probes),
#                nA (major allele CN), nB (minor allele CN), cn (total CN),
#                Ploidy, Aberrant Cell Fraction.
sequenza <- readRDS("Input/sequenza.rds")

# crc_cancer_genes : reference list of colorectal-cancer driver genes from NCG.
# Loaded here for downstream use; not plotted in this script.
load("Input/NCG_CRC_cancer_genes.Rdata")

# --- Export segment table ----------------------------------------------------
# Write the full Sequenza segment table to CSV for inspection / downstream use.
write.csv(sequenza, file = "Results/sequenza.results.csv", row.names = FALSE)

# --- Genome-wide CN plot -----------------------------------------------------
# Open a high-resolution PNG device (2000x1000 px, 200 dpi).
png(file = "Results/sequenza.png", w = 2000, h = 1000, res = 200)

y_limit    <- 5                        # maximum CN value shown on y-axis
len        <- sum(sequenza$nProbes)    # total probes across the genome (x-axis span)
twoColours <- TRUE                     # when TRUE, dim segments that exceed y_limit

# Set margins and character expansion factors for the base-R plot
par(mar = c(0.5, 5, 5, 0.5), cex = 0.4, cex.main = 3, cex.axis = 2.5)
ticks <- seq(0, y_limit, 1)           # y-axis tick marks at every integer CN

# Build run-length-encoding-style lists so we can iterate segment-by-segment.
# "lengths" = number of probes in each segment; "values" = CN of that segment.
A_rle <- list(lengths = sequenza$nProbes, values = sequenza$nA)   # major allele
B_rle <- list(lengths = sequenza$nProbes, values = sequenza$nB)   # minor allele

# Create an empty plot frame spanning the whole genome
plot(c(1, len), c(0, y_limit), type = "n", xaxt = "n", yaxt = "n",
     main = NULL, xlab = "", ylab = "")
axis(side = 2, at = ticks)                # draw y-axis CN labels
abline(h = ticks, col = "lightgrey", lty = 1)  # horizontal grid lines

colourMinor <- "yellow"   # colour for minor allele (B) segments
colourTotal <- "blue"     # colour for major allele (A) segments

# --- Draw minor-allele CN (yellow) rectangles --------------------------------
# Each segment is drawn as a thin horizontal rectangle centred on its CN value.
start <- 0
for (i in seq_along(B_rle$values)) {
  val  <- B_rle$values[i]    # minor-allele CN for this segment
  size <- B_rle$lengths[i]   # width in probes
  # Dim the colour if the CN exceeds the y-axis limit (avoids visual clutter)
  col_fill <- ifelse(twoColours & val >= y_limit,
                     adjustcolor(colourMinor, red.f = 0.75, green.f = 0.75, blue.f = 0.75),
                     colourMinor)
  rect(start, val - 0.07, start + size - 1, val + 0.07,
       col = col_fill, border = col_fill)
  start <- start + size
}

# --- Draw major-allele CN (blue) rectangles ----------------------------------
start <- 0
for (i in seq_along(A_rle$values)) {
  val  <- A_rle$values[i]    # major-allele CN for this segment
  size <- A_rle$lengths[i]
  col_fill <- ifelse(twoColours & val >= y_limit,
                     adjustcolor(colourTotal, red.f = 0.75, green.f = 0.75, blue.f = 0.75),
                     colourTotal)
  rect(start, val - 0.07, start + size - 1, val + 0.07,
       col = col_fill, border = col_fill)
  start <- start + size
}

# --- Chromosome boundary lines and labels ------------------------------------
# Compute cumulative probe count per chromosome to get boundary positions
chr.segs    <- ddply(sequenza, .(Chr), summarise, n = sum(nProbes))
chr.segs$cs <- cumsum(chr.segs$n)   # cumulative sum = right boundary of each chr

chrk_tot_len <- 0
abline(v = 0, lty = 1, col = "lightgrey")   # left genome boundary

for (i in 1:nrow(chr.segs)) {
  # Draw a vertical line at each chromosome boundary
  abline(v = chr.segs$cs[i], lty = 1, col = "lightgrey")

  # Place the chromosome label at the midpoint of each chromosome span
  chrk_tot_len_prev <- chrk_tot_len
  chrk_tot_len      <- chr.segs$cs[i]
  tpos <- (chrk_tot_len + chrk_tot_len_prev) / 2
  # Label chromosomes 1-22 numerically, then X and Y
  text(tpos, y_limit,
       ifelse(i < 23, sprintf("%d", i), ifelse(i == 23, "X", "Y")),
       pos = 1, cex = 2)
}

dev.off()   # close the PNG device and flush to disk
