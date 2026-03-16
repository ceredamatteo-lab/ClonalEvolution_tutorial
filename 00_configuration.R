## ============================================================================
## 00_configuration.R
## --------------------------------------------------------------------------
## Purpose : Install all R packages required by the Clonal Evolution tutorial.
##
## Packages fall into three categories:
##   1. CRAN      – general-purpose data wrangling and visualisation
##   2. Bioconductor – interval arithmetic (IRanges)
##   3. GitHub    – VAF-based mutation clustering (sciClone + its dependency bmm)
##
## Run this script once before executing the analysis pipeline (01–04).
## ============================================================================

# --- 1. CRAN packages --------------------------------------------------------
# RColorBrewer : colour palettes for plots
# plyr         : split-apply-combine operations (ddply)
# ggplot2      : grammar-of-graphics plotting
# reshape2     : wide <-> long data reshaping (melt / dcast)
# ggrepel      : non-overlapping text labels on ggplot2 plots
# devtools     : install packages from GitHub repositories
pkgs <- c('RColorBrewer', 'plyr', 'ggplot2', 'reshape2', 'ggrepel', 'devtools')

# Check which packages are already installed and install only the missing ones
lib <- installed.packages()
installed <- pkgs %in% rownames(lib)
not_installed <- which(!installed)
if (length(not_installed) > 0) {
  for (i in not_installed) install.packages(pkgs[i])
}

# --- 2. Bioconductor packages ------------------------------------------------
# BiocManager replaces the deprecated biocLite() installer (Bioconductor >= 3.8).
# IRanges provides efficient interval (range) operations used internally.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("IRanges")

# --- 3. GitHub packages -------------------------------------------------------
# bmm      : Bayesian mixture model engine (dependency of sciClone)
# sciClone : clusters somatic variants into subclones using VAF distributions
#            and copy-number-aware filtering (Miller et al., PLoS Comput Biol 2014)
library(devtools)
install_github("genome/bmm")
install_github("genome/sciClone")
