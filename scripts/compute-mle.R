# compute-mle.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Jan 8, 2020
#
# Run function to compute poisson lognormal mle parameters

library(here)
source(here("scripts", "quminorm.R"))

# Load filtered counts
filtered_file_name = here("salmon_quants", "dat_filtered.rds")
if(file.exists(filtered_file_name)){
  dat_filtered = readRDS(filtered_file_name)
}

# Find mle parameters mu and sig
if(use_matrix){
  m = round(dat_filtered)
  res = poilog_mle_matrix(m)
  saveRDS(res, here("mle_results", "all.rds"))
} else {
  #$ -t 1-534
  #$ -tc 50
  i = as.integer(Sys.getenv("SGE_TASK_ID")) 
  x = dat_filtered[, i]
  x = round(x) # counts must be integers
  res = poilog_mle(x)
  saveRDS(res, here("mle_results", paste0(i, ".rds")))
}