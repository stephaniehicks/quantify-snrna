# Reproducible workflows for the analysis of snRNA-seq data


## Project: Mouse Cortex 

The `/mouse_cortex/` folder contains code to download and process data from this [paper](https://www.nature.com/articles/s41587-020-0465-8). Data processing and analysis are described below (more details in `/mouse_cortex/code/workflow.Rmd`).

### Download data

#### Download GEO metadata

First, you need to manually download the `SraRunInfo.csv` file. 

* Go here: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA545730. 
* Then click on `SRA Experiments`. Click `Send to`. Choose `File`. 
* Change Format to `RunInfo`. Click `Create File`. 
* This will download a file called `SraRunInfo.csv`.

Then download and create GEO metadata with the shell script `download-geo-data.sh`. 

#### Download SRA files

Install the SRA toolkit. Here are some instructions on how to do so:

1. https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/
2. https://github.com/ncbi/sra-tools/wiki/05.-Toolkit-Configuration (old version is https://ncbi.github.io/sra-tools/install_config.html)

Notes:
* You will need to add the Toolkit functions to your PATH variable (https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit). This can be done by adding ~/.../sratoolkit.2.10.6-centos_linux64/bin to the PATH variable in ~/.bash_profile. Check that it worked with `which fastq-dump`.
* You will also want to specify the default location to download the SRA files to using the toolkit.

Once you have installed the toolkit, run the shell script `download-sra-data.sh`. An example of a description of an SRA file is https://www.ncbi.nlm.nih.gov/sra/SRX5943587. Note that the file name Cortex1.CCJ15ANXX.10X_2B indicates the cortex (Cortex1), flow cell (CCJ15ANXX), sequencing method (10X), and the lane (2B).

#### Extract FASTQ files

References:
1. https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump
2. https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump

After downloading the SRA files, convert them to FASTQ files using the shell script `extract-fastq.sh`.

#### Download GENCODE files

Run the shell script `download-gencode-files.sh`. This downloads GENCODE files and creates the files necessary for *alevin*. Reference files were obtained for [mouse release M25](https://www.gencodegenes.org/mouse/release_M25.html). 
Note that the decoy-aware transcriptomes are created directly in the shell script `download-gencode-files.sh` and not the R script it calls. Details about these steps are in `quantify-salmon.Rmd`.

### Run alevin

Build the reference index using `build-index-salmon.sh`. This process takes a few hours for each pipeline. Then run alevin using `run-alevin.sh`.

For now, you need to manually specify which pipeline you are running in both of these scripts.

### Clean data

Read in the quant files from alevin into R using `run-tximeta.sh`. Since alevin import only supports one file at a time currently, this creates a separate SummarizedExperiment object for each library. 

Combine the 2 SummarizedExperiment objects together using `combine-se.sh`. 

Then run quality control, PCA, and add cell type labels with `save-sce.R`. This also converts the SummarizedExperiment object into a SingleCellExperiment object.

### Exploratory analysis and plots

Exploratory work and plots is done in `snuc-analysis-bioc.Rmd`. Publication figures are generated using `distribution-plots.sh`, `libsize-plots.sh`, `celltype-plots.sh`, and `pca-plots.sh`, and tweaked with `generate-figures.Rmd`. 

## Project: snRNA-seq from pediatric HGG tumors (Filbin)

### files 
1. `/human_cortex/code/quantify-salmon.Rmd` contains the code to:
- download files needed for Salmon
    - save a combined pre-mRNA and mRNA fasta and gtf file 
    - install salmon v.1.0.0
    - build the salmon index with pre-mRNA and mRNA (see `/human_cortex/code/build-index-salmon.sh`)
    - run salmon for quantification of counts for all tumors (see `/human_cortex/code/run-salmon.sh`)
    
2. `/human_cortex/code/snuc-analysis-bioc.Rmd` contains the code to: 
- create `SummarizedExperiment` object with the `tximeta` R/Bioconductor package
    - convert to a `SingleCellExperiment` object
    - quality control and preprocessing using `scater`

## Authors

* Stephanie Hicks (shicks19@jhu.edu)
* Albert Kuo (albertkuo@jhu.edu)