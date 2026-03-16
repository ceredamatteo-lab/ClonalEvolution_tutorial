## ============================================================================
## 03_clone_composition.R
## --------------------------------------------------------------------------
## Purpose : Classify each somatic variant as monoclonal, biclonal, or
##           polyclonal based on its clonality score, then visualise the
##           composition as a stacked horizontal bar chart per patient.
##
## Thresholds (following Caravagna et al.):
##   - Monoclonal  : clonality >= 0.80  (present in virtually all cancer cells)
##   - Biclonal    : 0.35 <= clonality < 0.80  (present in a major sub-population)
##   - Polyclonal  : clonality < 0.35  (present in a minor sub-population)
##
## Input files:
##   - Results/SNVS.clonality.Rdata  (output of 02_measure_clonality.R)
##
## Output files:
##   - Results/CloneComposition.pdf  (stacked bar chart of clone fractions)
## ============================================================================

library(plyr)          # ddply for per-patient summarisation
library(reshape2)      # melt: wide -> long format for ggplot stacking
library(ggplot2)       # bar chart
library(ggrepel)       # (loaded for consistency; not used directly here)
library(RColorBrewer)  # (loaded for consistency; custom palette defined below)

# --- Load clonality-annotated variants produced by step 02 -------------------
load("Results/SNVS.clonality.Rdata")   # provides: snvs (with frq.tc, clonality)

# --- Colour palette for the three clone categories ---------------------------
# M = monoclonal (blue), B = biclonal (pink), P = polyclonal (yellow)
color_clone_composition <- c(
  'M' = rgb(0,   162, 205, maxColorValue = 255),
  'B' = rgb(243, 130, 153, maxColorValue = 255),
  'P' = rgb(223, 207,   0, maxColorValue = 255)
)

# --- Minimal ggplot2 theme (clean, no gridlines) -----------------------------
theme_cloneR <- function(base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.background = element_blank(),
          panel.border     = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}

# --- Define clonality thresholds ---------------------------------------------
upper <- 0.80   # boundary between monoclonal and biclonal
lower <- 0.35   # boundary between biclonal and polyclonal

# Store thresholds in the data.frame (used inside ddply via unique())
snvs$upper <- upper
snvs$lower <- lower

# --- Per-patient summary: count and fraction in each category ----------------
y <- ddply(snvs, .(Patient), summarise,
           n            = length(clonality),                                   # total variants
           n_monoclonal = sum(clonality >= unique(upper)),                     # clonality >= 0.80
           n_biclonal   = sum(clonality <  unique(upper) & clonality >= unique(lower)),  # 0.35-0.80
           n_polyclonal = sum(clonality <  unique(lower)))                     # clonality <  0.35

# Convert counts to proportions (0-1) for stacking
y$monoclonal <- y$n_monoclonal / y$n
y$biclonal   <- y$n_biclonal   / y$n
y$polyclonal <- y$n_polyclonal / y$n

# Assign a single-letter dominant-composition code based on the largest fraction
code <- c("M", "B", "P")
names(code) <- c('monoclonal', 'biclonal', 'polyclonal')
y$composition <- code[names(which.max(y[, c('monoclonal', 'biclonal', 'polyclonal')]))]

# --- Plotting function -------------------------------------------------------
# Produces a horizontal stacked bar chart showing the proportion of variants
# in each clonality category for every patient.
clone.composition.plot <- function(x, cl = color_clone_composition) {

  # Map colour names to the full category labels expected by scale_fill_manual
  names(cl) <- c("monoclonal", "biclonal", "polyclonal")

  # Ensure composition is an ordered factor for consistent legend order
  x$composition <- factor(x$composition, levels = c("M", "B", "P"))

  if (!is.null(x)) {
    # Build a patient -> composition lookup (useful for facet labelling)
    mapper <- as.list(x$composition)
    names(mapper) <- x$Patient
    map_labeller <- function(variable, value) {
      return(mapper[value])
    }

    # Reshape from wide (one column per category) to long format for geom_bar
    m <- melt(x[, c('Patient', 'polyclonal', 'biclonal', 'monoclonal')],
              id.vars = c("Patient"))
    colnames(m)[2] <- 'composition'

    # Build the stacked bar plot
    p <- ggplot(m, aes(x = Patient, y = value, fill = composition)) +
      geom_bar(width = 0.5, stat = "identity") +
      # Draw a solid left-side axis line from 0 to 1
      geom_segment(aes(x = -Inf, xend = -Inf, y = 0, yend = 1), col = "black") +
      ylab("Alterations (%)") + xlab("") +
      scale_fill_manual(
        values = cl,
        guide  = guide_legend(title = NULL),
        labels = c("Clonality<35%", "35%<Clonality<80%", "Clonality>80%")
      ) +
      scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
      theme_cloneR() +
      theme(panel.background = element_blank(),
            legend.position  = "top",
            legend.key       = element_rect(size = 1.5, color = 'white'),
            legend.text      = element_text(color = "black", size = 10),
            axis.text        = element_text(color = "black", size = 10),
            axis.text.y      = element_blank(),
            axis.ticks.y     = element_blank(),
            strip.background = element_rect(fill = NA, colour = NA),
            strip.text.y     = element_text(angle = 0, size = 12, colour = "black")) +
      # NOTE: In ggplot2 only the last coord_* layer takes effect, so
      # coord_equal() below is overridden by coord_flip().  Both are kept
      # to preserve the original pipeline behaviour.
      coord_equal(1 / 0.1) +
      coord_flip() +
      # Re-draw bars with a black outline on top of the filled bars
      geom_bar(width = 0.5, stat = "identity", color = "black",
               show.legend = FALSE)

    return(p)
  } else {
    return(NULL)
  }
}

# --- Generate and save the composition plot ----------------------------------
pdf(file = "Results/CloneComposition.pdf", h = 2, w = 6)
clone.composition.plot(y)
dev.off()
