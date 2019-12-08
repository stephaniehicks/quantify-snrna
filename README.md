# Reproducible workflows for the analysis of snRNA-seq data

## Summary of relevant files

1. `/scripts/quantify-salmon.Rmd` contains the code to:

    - build the salmon index with pre-mRNA and mRNA (`/scripts/build-index-salmon.sh`)
    - run salmon for quantification of counts for all tumors (`/scripts/run-salmon.sh`)
    - create `SummarizedExperiment` object

2. 



## Authors

* Stephanie Hicks (shicks19@jhu.edu)
* Albert Kuo (albertkuo@jhu.edu)