# pca-plots.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Sep 24, 2020
#
# Make all PCA plots


suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(scran)
  library(scater)
  library(BiocSingular)
  library(SingleCellExperiment)
})


# Read in SingleCellExperiment objects
run_number = "all" # give run_number or "all" for all of them together
sce_ls = list()
sce_ls[["transcripts"]] = readRDS(here("mouse_cortex", "salmon_quants", "transcripts_pipeline", paste0("sce_", run_number, ".rds")))
sce_ls[["preandmrna"]] = readRDS(here("mouse_cortex", "salmon_quants", "preandmrna_pipeline", paste0("sce_", run_number, ".rds")))
sce_ls[["introncollapse"]] = readRDS(here("mouse_cortex", "salmon_quants", "introncollapse_pipeline", paste0("sce_", run_number, ".rds")))
sce_ls[["intronseparate"]] = readRDS(here("mouse_cortex", "salmon_quants", "intronseparate_pipeline", paste0("sce_", run_number, ".rds")))


#############################
# Plot PCA with Ding labels #
#############################
p = plotReducedDim(sce_ls[["transcripts"]], "PCA", colour_by = "ding_labels")
ggsave(file = here(paste0("./mouse_cortex/plots/pcading_plot_transcripts.png")), plot = p)
saveRDS(p, here(paste0("./mouse_cortex/plots/pcading_plot_transcripts.rds")))

p = plotReducedDim(sce_ls[["preandmrna"]], "PCA", colour_by = "ding_labels")
ggsave(file = here(paste0("./mouse_cortex/plots/pcading_plot_preandmrna.png")), plot = p)
saveRDS(p, here(paste0("./mouse_cortex/plots/pcading_plot_preandmrna.rds")))

p = plotReducedDim(sce_ls[["introncollapse"]], "PCA", colour_by = "ding_labels")
ggsave(file = here(paste0("./mouse_cortex/plots/pcading_plot_introncollapse.png")), plot = p)
saveRDS(p, here(paste0("./mouse_cortex/plots/pcading_plot_introncollapse.rds")))

p = plotReducedDim(sce_ls[["intronseparate"]], "PCA", colour_by = "ding_labels")
ggsave(file = here(paste0("./mouse_cortex/plots/pcading_plot_intronseparate.png")), plot = p)
saveRDS(p, here(paste0("./mouse_cortex/plots/pcading_plot_intronseparate.rds")))