# download-gencode-files.R
# -----------------------------------------------------------------------------
# Author:             Albert Kuo
# Date last modified: Oct 8, 2020
#
# Download GENCODE files and generate files for the pipelines

library(here)

# Create directories
if(!dir.exists(here("mouse_cortex", "salmon_files"))){
  dir.create(here("salmon_files", "human"))
}
if(!dir.exists(here("mouse_cortex", "salmon_quants"))){
  dir.create(here("mouse_cortex", "salmon_quants"))
}

# Download GENCODE Files
# download GENCODE primary assembly fasta file
if(!file.exists(here("mouse_cortex", "salmon_files", "GRCm38.primary_assembly.genome.fa.gz"))){
  tar_gz_file <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz"
  download.file(tar_gz_file, 
                destfile = here("mouse_cortex", "salmon_files", "GRCm38.primary_assembly.genome.fa.gz"), 
                method = "wget")
}

# download GENCODE transcripts fasta file (deprecated)
# if(!file.exists(here("mouse_cortex", "salmon_files", "gencode.vM25.transcripts.fa.gz"))){
#   tar_gz_file <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.transcripts.fa.gz"
#   download.file(tar_gz_file,
#                 destfile = here("mouse_cortex", "salmon_files", "gencode.vM25.transcripts.fa.gz"),
#                 method = "wget")
# }

# download GENCODE gtf file
if(!file.exists(here("mouse_cortex", "salmon_files", "gencode.vM25.annotation.gtf.gz"))){
  tar_gz_file <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"
  download.file(tar_gz_file, 
                destfile = here("mouse_cortex", "salmon_files", "gencode.vM25.annotation.gtf.gz"), 
                method = "wget")
}


# Make files for pipelines (code adapted from https://github.com/csoneson/rna_velocity_quant)
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(Biostrings)
  library(rtracklayer)
  library(GenomicFeatures)
  library(BSgenome)
})
source(here("mouse_cortex", "code", "quantify-salmon-helpers.R")) 


#######################
# transcripts pipline #
#######################
## FASTA file
# Gtf path
gtf_file <- here("mouse_cortex", "salmon_files", "gencode.vM25.annotation.gtf.gz")

# Read genomic (DNA) sequence from FASTA file
genome_fasta <- here("mouse_cortex", "salmon_files", "GRCm38.primary_assembly.genome.fa.gz") 
genome <- Biostrings::readDNAStringSet(genome_fasta)
names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1) # creates chr1, etc

# Extract transcript (tx) sequences (takes a few minutes)
tx <- extractTxSeqs(gtf = gtf_file, genome = genome, type = "spliced")

# Write FASTA file
Biostrings::writeXStringSet(tx, file = here("mouse_cortex", "salmon_files", "gencode.vM25.transcripts.mouse.fa.gz"), compress = T) # Compresses fa while saving, rather than compress later (takes a few min)

## tx2gene
# Read gtf
gtf <- rtracklayer::import(here("mouse_cortex", "salmon_files",
                                "gencode.vM25.annotation.gtf.gz")) 
gtftx <- subset(gtf, type == "transcript")
gtfex <- subset(gtf, type == "exon")

df <- data.frame(gtftx, stringsAsFactors = FALSE) %>%
  dplyr::select(transcript_id, seqnames, start, end, strand, source, 
                gene_id, gene_type, gene_name, level, havana_gene, transcript_type,
                transcript_name, transcript_support_level, tag, havana_transcript) %>%
  dplyr::left_join(data.frame(gtfex, stringsAsFactors = FALSE) %>%
                     dplyr::group_by(transcript_id) %>% 
                     dplyr::summarize(transcript_length = sum(width)),
                   by = "transcript_id")

# Write table as txt and rds
write.table(df %>% dplyr::select(transcript_id, gene_id), 
            file = here("mouse_cortex", "salmon_files", "gencode.vM25.transcripts.tx2gene.mouse.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE, 
            col.names = FALSE)
saveRDS(df, file = here("mouse_cortex", "salmon_files", "gencode.vM25.transcripts.tx2gene.mouse.rds"))


#######################
# preandmrna pipline ##
#######################
## FASTA file
# Gtf path
gtf_file <- here("mouse_cortex", "salmon_files", "gencode.vM25.annotation.gtf.gz")

# Read genomic (DNA) sequence from FASTA file
genome_fasta <- here("mouse_cortex", "salmon_files", "GRCm38.primary_assembly.genome.fa.gz") 
genome <- Biostrings::readDNAStringSet(genome_fasta)
names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1) # creates chr1, etc

# Extract transcript (tx) and pre-mRNA (premrna) sequences (takes a few minutes)
tx <- extractTxSeqs(gtf = gtf_file, genome = genome, type = "spliced")
premrna <- extractTxSeqs(gtf = gtf_file, genome = genome, type = "unspliced")
names(premrna) <- paste0(names(premrna), "_unspliced")

# Combine mRNA and pre-mRNA sequences
preandmrna <- c(tx, premrna)

# Write FASTA file
Biostrings::writeXStringSet(preandmrna, file = here("mouse_cortex", "salmon_files", "gencode.vM25.preandmrna.mouse.fa.gz"), compress = T) # Compresses fa while saving, rather than compress later (takes a while ~30 min)


## tx2gene
# Generate tx2gene table for mRNA and pre-mRNA transcripts
t2g <- readRDS(here("mouse_cortex", "salmon_files", "gencode.vM25.transcripts.tx2gene.mouse.rds"))
t2gpre <- t2g %>% dplyr::mutate(transcript_id = paste0(transcript_id, "_unspliced"))
t2g <- rbind(t2g, t2gpre)

# Write table as txt
write.table(t2g %>% dplyr::select(transcript_id, gene_id), 
            file = here("mouse_cortex", "salmon_files", "gencode.vM25.preandmrna.tx2gene.mouse.txt"), 
            row.names = FALSE, col.names = FALSE, 
            sep = "\t", quote = FALSE)


###########################
# introncollapse pipline ##
###########################
## FASTSA file
# Gtf path
gtf_file <- here("mouse_cortex", "salmon_files", "gencode.vM25.annotation.gtf.gz")

# Read genomic (DNA) sequence from FASTA file
genome_fasta <- here("mouse_cortex", "salmon_files", "GRCm38.primary_assembly.genome.fa.gz") 
genome <- Biostrings::readDNAStringSet(genome_fasta)
names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1) # creates chr1, etc

# Extract transcript (tx) and intron sequences (takes a few minutes)
tx <- extractTxSeqs(gtf = gtf_file, genome = genome, type = "spliced")
isoform_action = "collapse"
flanklength = 50 # read length
intr <- extractIntronSeqs(gtf = gtf_file, genome = genome, type = isoform_action, 
                          flanklength = flanklength,
                          joinOverlappingIntrons = FALSE)

# Combine mRNA and pre-mRNA sequences
intronandmrna <- c(tx, intr)

# Write FASTA file
Biostrings::writeXStringSet(intronandmrna, file = here("mouse_cortex", "salmon_files",
                                                       paste0("gencode.vM25.intron", isoform_action, ".mouse.fa.gz")), compress = TRUE)

## tx2gene
# Generate transcript/intron-to-gene mapping
gtfdf <- as.data.frame(rtracklayer::import(here("mouse_cortex", "salmon_files",
                                                "gencode.vM25.annotation.gtf.gz")))
t2gtx <- gtfdf %>% dplyr::filter(type == "transcript") %>%
  dplyr::select(transcript_id, gene_id) %>%
  dplyr::distinct()

if (isoform_action == "collapse") {
  ## Intron names already contain gene name
  t2gin <- data.frame(intr = names(intr),
                      gene = gsub("\\-I[0-9]*$", "", names(intr)),
                      stringsAsFactors = FALSE)
} else if (isoform_action == "separate") {
  ## Intron names contain transcript name
  t2gin <- data.frame(intr = names(intr),
                      transcript_id = gsub("\\-I[0-9]*$", "", names(intr)),
                      stringsAsFactors = FALSE) %>%
    dplyr::left_join(t2gtx, by = "transcript_id") %>%
    dplyr::select(intr, gene_id)
} else {
  stop("Unknown isoform_action")
}
colnames(t2gin) <- colnames(t2gtx)
t2g <- rbind(t2gtx, t2gin)

# Write table as txt
write.table(t2g, file = here("mouse_cortex", "salmon_files", paste0("gencode.vM25.intron", isoform_action, ".tx2gene.mouse.txt")), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


###########################
# intronseparate pipline ##
###########################
## FASTSA file
# Gtf path
gtf_file <- here("mouse_cortex", "salmon_files", "gencode.vM25.annotation.gtf.gz")

# Read genomic (DNA) sequence from FASTA file
genome_fasta <- here("mouse_cortex", "salmon_files", "GRCm38.primary_assembly.genome.fa.gz") 
genome <- Biostrings::readDNAStringSet(genome_fasta)
names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1) # creates chr1, etc

# Extract transcript (tx) and intron sequences (takes a few minutes)
tx <- extractTxSeqs(gtf = gtf_file, genome = genome, type = "spliced")
isoform_action = "separate"
flanklength = 50 # read length
intr <- extractIntronSeqs(gtf = gtf_file, genome = genome, type = isoform_action, 
                          flanklength = flanklength,
                          joinOverlappingIntrons = FALSE)

# Combine mRNA and pre-mRNA sequences
intronandmrna <- c(tx, intr)

# Write FASTA file
Biostrings::writeXStringSet(intronandmrna, file = here("mouse_cortex", "salmon_files",
                                                       paste0("gencode.vM25.intron", isoform_action, ".mouse.fa.gz")), compress = TRUE)

## tx2gene
# Generate transcript/intron-to-gene mapping
gtfdf <- as.data.frame(rtracklayer::import(here("mouse_cortex", "salmon_files",
                                                "gencode.vM25.annotation.gtf.gz")))
t2gtx <- gtfdf %>% dplyr::filter(type == "transcript") %>%
  dplyr::select(transcript_id, gene_id) %>%
  dplyr::distinct()

if (isoform_action == "collapse") {
  ## Intron names already contain gene name
  t2gin <- data.frame(intr = names(intr),
                      gene = gsub("\\-I[0-9]*$", "", names(intr)),
                      stringsAsFactors = FALSE)
} else if (isoform_action == "separate") {
  ## Intron names contain transcript name
  t2gin <- data.frame(intr = names(intr),
                      transcript_id = gsub("\\-I[0-9]*$", "", names(intr)),
                      stringsAsFactors = FALSE) %>%
    dplyr::left_join(t2gtx, by = "transcript_id") %>%
    dplyr::select(intr, gene_id)
} else {
  stop("Unknown isoform_action")
}
colnames(t2gin) <- colnames(t2gtx)
t2g <- rbind(t2gtx, t2gin)

# Write table as txt
write.table(t2g, file = here("mouse_cortex", "salmon_files", paste0("gencode.vM25.intron", isoform_action, ".tx2gene.mouse.txt")), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
