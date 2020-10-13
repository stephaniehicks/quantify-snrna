# save-sce.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Oct 13, 2020
#
# Quality control, run PCA, add cell type labels, and save as SingleCellExperiment 

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

# Combine counts for same cell barcode
for(i in seq_along(se_ls)){
  # Sum up counts
  m = assay(se, "counts")
  colnames(m) = paste0(colData(se)$cortex, "_", colnames(m)) # colnames format is Cortex1_CGACTTCAGTCTCGGC
  m = colsum(m, colnames(m))                                 # Sum up counts for same column name
  
  # Create new colData
  col_data = colData(se) %>% 
    as_tibble() %>%
    select(-Run, -flow_cell, -lane) %>%
    mutate(cell_barcode = paste0(cortex, "_", rownames(colData(se)))) %>%
    distinct() %>%
    left_join(data.frame(cell_barcode = colnames(m)), . , by = "cell_barcode") # Arrange in same order as m
  
  # Create new SE
  se_new = SummarizedExperiment(assays = SimpleList(counts = m),
                                colData = col_data,
                                rowData = rowData(se))
  
  se_ls[[i]] = se_new
}

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

# Cell type labels (Ding)
meta_ding = read_tsv(here("mouse_cortex", "files", "meta_combined.txt"))
meta_ding_10x = meta_ding %>%
  filter(Method == "10x Chromium") %>%
  mutate(cell_barcode = paste0(Experiment, "_", gsub("Cortex.*10xChromium", "", NAME)))

for(i in seq_along(sce_ls)){
  sce = sce_ls[[i]]
  
  match_cols = match(colnames(sce), meta_ding_10x$cell_barcode)
  colData(sce)$ding_labels = meta_ding_10x$CellType[match_cols]
  
  sce_ls[[i]] = sce
}

# Cell type labels (singleR)
# Built-in reference mouse dataset
library(SingleR)
mouse_se = MouseRNAseqData() # inclusion in singleR is deprecated, obtain from celldex

for(i in seq_along(sce_ls)){
  sce = sce_ls[[i]]
  
  # Convert gene names
  test_counts = assay(sce, "logcounts")
  tx2gene = readRDS(here("mouse_cortex", "salmon_files", "gencode.vM25.transcripts.tx2gene.mouse.rds"))
  match_rows = match(rownames(test_counts), tx2gene$gene_id)
  rownames(test_counts) = tx2gene$gene_name[match_rows]
  
  pred_celltypes = SingleR(test = test_counts, 
                           ref = mouse_se,
                           labels = mouse_se$label.fine)
  
  all.markers <- metadata(pred_celltypes)$de.genes
  
  # Save SingleR output
  saveRDS(pred_celltypes, here("mouse_cortex", "salmon_quants", paste0(names(sce_ls)[[i]], "_pipeline/singler_results.rds")))
  
  # Add to colData
  colData(sce)$singleR_labels = pred_celltypes$labels
  colData(sce)$singleR_labels_pruned = pred_celltypes$pruned.labels
  colData(sce)$singleR_main_labels = colData(sce) %>% 
    as_tibble() %>%
    select(singleR_labels) %>%
    left_join(., as_tibble(colData(mouse_se)) %>% distinct(), by = c("singleR_labels" = "label.fine")) %>%
    select(label.main)
  
  sce_ls[[i]] = sce
}


# Save SCE with PCA
saveRDS(sce_ls[["transcripts"]], here("mouse_cortex", "salmon_quants", "transcripts_pipeline", paste0("sce_", run_number, ".rds")))
saveRDS(sce_ls[["preandmrna"]], here("mouse_cortex", "salmon_quants", "preandmrna_pipeline", paste0("sce_", run_number, ".rds")))
saveRDS(sce_ls[["introncollapse"]], here("mouse_cortex", "salmon_quants", "introncollapse_pipeline", paste0("sce_", run_number, ".rds")))
saveRDS(sce_ls[["intronseparate"]], here("mouse_cortex", "salmon_quants", "intronseparate_pipeline", paste0("sce_", run_number, ".rds")))
