# save-sce.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Sep 9, 2020
#
# Quality control, run PCA and save as SingleCellExperiment 

suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(tictoc)
  library(scran)
  library(scater)
  library(BiocSingular)
  library(SingleCellExperiment)
})

run_number = "all" # give run_number or "all" for all of them together
se_ls = list()
se_ls[["transcripts"]] = readRDS(here("mouse_cortex", "salmon_quants", "transcripts_pipeline", paste0("se_", run_number, ".rds")))
se_ls[["preandmrna"]] = readRDS(here("mouse_cortex", "salmon_quants", "preandmrna_pipeline", paste0("se_", run_number, ".rds")))
se_ls[["introncollapse"]] = readRDS(here("mouse_cortex", "salmon_quants", "introncollapse_pipeline", paste0("se_", run_number, ".rds")))
se_ls[["intronseparate"]] = readRDS(here("mouse_cortex", "salmon_quants", "intronseparate_pipeline", paste0("se_", run_number, ".rds")))

# Quality control
sce_ls = list()
for(i in seq_along(se_ls)){
  se = se_ls[[i]]
  # Remove empty rows
  non_empty_rows <- which(rowSums(assay(se, "counts")) > 1)
  se_subset <- se[non_empty_rows, ]
  
  # Remove low-quality cells
  df = perCellQCMetrics(se_subset)
  reasons = quickPerCellQC(df) # DataFrame of logical values
  colSums(as.matrix(reasons))
  
  # Convert to SingleCellExperiment
  sce = as(se_subset[, !reasons$discard], "SingleCellExperiment")
  colData(sce)$lib_size = colSums(counts(sce))
  sce_ls[[names(se_ls)[i]]] = sce
}
rm(se_ls, se, se_subset)

# PCA
for(i in seq_along(sce_ls)){
  sce = sce_ls[[i]]
  # Normalize and take log
  clust = quickCluster(sce)
  sce = computeSumFactors(sce, clusters = clust)
  sce = logNormCounts(sce)
  
  set.seed(1)
  tic("approx PCA")
  sce = runPCA(sce, exprs_values = "logcounts", ntop = ncol(sce), BSPARAM = RandomParam())
  toc()
  
  sce_ls[[i]] = sce
}

# Save SCE with PCA
saveRDS(sce_ls[["transcripts"]], here("mouse_cortex", "salmon_quants", "transcripts_pipeline", paste0("sce_", run_number, ".rds")))
saveRDS(sce_ls[["preandmrna"]], here("mouse_cortex", "salmon_quants", "preandmrna_pipeline", paste0("sce_", run_number, ".rds")))
saveRDS(sce_ls[["introncollapse"]], here("mouse_cortex", "salmon_quants", "introncollapse_pipeline", paste0("sce_", run_number, ".rds")))
saveRDS(sce_ls[["intronseparate"]], here("mouse_cortex", "salmon_quants", "intronseparate_pipeline", paste0("sce_", run_number, ".rds")))
