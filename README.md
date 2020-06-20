# Reproducible workflows for the analysis of snRNA-seq data


## Project: Mouse Cortex 

The `/mouse_cortex/` folder contains code to download and process data from this [paper](https://www.nature.com/articles/s41587-020-0465-8): 

### files 

    - `/mouse_cortex/download-geo-data.Rmd`: Download the SRA files from GEO
        - `/mouse_cortex/download-geo-data.sh`: shell script to prefetch SRA files
        - `/mouse_cortex/extract-fastq.sh`: shell script to convert SRA to fastq
    - `/mouse_cortex/quantify-salmon.Rmd`: quantification with alevin from salmon
        - `/mouse_cortex/build-index-salmon.sh`: shell script to build salmon index
        - `/mouse_cortex/run-alevin.sh`: shell script to run alevin

### reference generation

All reference files were obtained from GENCODE, [mouse release M25](https://www.gencodegenes.org/mouse/release_M25.html). 
The following indices were prepared:

#### Salmon index (using decoys)

```
salmon index -t gentrome_transcripts.fa.gz \
                -d decoys.txt \
                -i salmon_transcripts_index --gencode --threads 4
```


## Project: snRNA-seq from pediatric HGG tumors (Filbin)

### files 
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