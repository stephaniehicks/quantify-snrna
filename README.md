# Reproducible workflows for the analysis of snRNA-seq data

## Summary of relevant files

1. `/scripts/quantify-salmon.Rmd` contains the code to:

    - download files needed for Salmon
    - save a combined pre-mRNA and mRNA fasta and gtf file 
    - install salmon v.1.0.0
    - build the salmon index with pre-mRNA and mRNA (see `/scripts/build-index-salmon.sh`)
    - run salmon for quantification of counts for all tumors (see `/scripts/run-salmon.sh`)

2. `/scripts/snuc-analysis-bioc.Rmd` contains the code to: 

    - create `SummarizedExperiment` object with the `tximeta` R/Bioconductor package
    - convert to a `SingleCellExperiment` object
    - quality control and preprocessing using `scater`

## Authors

* Stephanie Hicks (shicks19@jhu.edu)
* Albert Kuo (albertkuo@jhu.edu)