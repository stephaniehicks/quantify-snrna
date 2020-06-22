# compute-mle.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Feb 6, 2020
#
# Run function to compute poisson lognormal mle parameters

library(here)
library(SingleCellExperiment)
source(here("scripts", "quminorm.R"))

# Load filtered counts
dataset = "mrna"
filtered_file_name = here("salmon_quants", paste0("sce_filtered_", dataset, ".rds"))
sce_filtered = readRDS(filtered_file_name)
dat_filtered = counts(sce_filtered)

# Find mle parameters mu and sig
use_matrix = F
if(use_matrix){
  m = round(dat_filtered)
  res = poilog_mle_matrix(m)
  saveRDS(res, here(paste0("mle_results_", dataset), "all.rds"))
} else {
  #$ -t 1-576
  #$ -tc 50
  i = as.integer(Sys.getenv("SGE_TASK_ID")) 
  x = dat_filtered[, i]
  x = round(x) # counts must be integers
  res = poilog_mle(x)
  saveRDS(res, here(paste0("mle_results_", dataset), paste0(i, ".rds")))
}
