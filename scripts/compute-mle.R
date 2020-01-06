# compute-mle.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Jan 6, 2020
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
m = round(dat_filtered) # Must be integers
chunksize = 50          # chunksize = 50 cells is what Will used
res = parcolapply(m, poilog_mle_matrix, chunksize = chunksize)

saveRDS(res, here("salmon_quants", "mle_results.rds"))