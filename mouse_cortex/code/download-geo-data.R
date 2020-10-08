# download-geo-data.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Oct 8, 2020
#
# Download GEO data

library(here)
library(tidyverse)

# Download phenotype table
gse <- GEOquery::getGEO("GSE132044")
sapply(gse, dim) # Number of rows per entry
gse <- gse[[2]] # this contains all the mouse cortex data on 4 platforms
pdata <- Biobase::pData(gse) # Get phenotype table
pdata$Experiment <- sapply(stringr::str_split(pdata$relation.1, "="), tail, 1)

# To get SRR IDs, go to https://www.ncbi.nlm.nih.gov/bioproject/PRJNA545730. 
# Then click on `SRA Experiments`. Click `Send to`. Choose `File`. 
# Change Format to `RunInfo`. Click `Create File`. 
# This will download a file called `SraRunInfo.csv`. Read this file in. 
sra <- readr::read_csv(here("mouse_cortex", "files", "SraRunInfo.csv"))
pdata <- dplyr::left_join(pdata, sra, by = "Experiment")

# Save phenotype table
readr::write_csv(pdata, path = here("mouse_cortex", "files", "pdata.csv"))

# Shows there are four types of platforms used for library construction
table(pdata$characteristics_ch1.3)

# Select subset of phenotype table
sra_meta <- pdata %>% 
  select(Run, Experiment, title, geo_accession, type, source_name_ch1, organism_ch1, starts_with("characteristics"), molecule_ch1, taxid_ch1, description, description.1,
         platform_id, instrument_model, library_selection, library_source, library_strategy,
         `library construction:ch1`, `age:ch1`, `Sex:ch1`, `strain:ch1`, spots, bases, spots_with_mates, avgLength, size_MB, AssemblyName, SRAStudy, BioProject)

# All SRA files
write_csv(sra_meta, path = here("mouse_cortex", "files", "SRA_geo_metadata.csv"))

# 16 SRA files (10x)
pdata_10x <- pdata[pdata$characteristics_ch1.3 == 
                     "library construction: 10x Chromium (v2)",] 
write.table(pdata_10x$Run, file = here("mouse_cortex", "files", "SRR_files_10x.txt"), 
            quote= FALSE,row.names = FALSE, col.names = FALSE)
write.table(pdata_10x$download_path, file = here("mouse_cortex", "files", "SRR_paths_10x.txt"), 
            quote= FALSE,row.names = FALSE, col.names = FALSE)