# quantify-salmon-helpers.R
# -----------------------------------------------------------------------------
# Author:             Charlotte Soneson, Albert Kuo
# Date last modified: Jun 26, 2020
#
# Functions for quantification from https://github.com/csoneson/rna_velocity_quant/

#' Extract (spliced or unspliced) transcript sequences
#'
#' @param gtf The path to a gtf file
#' @param genome A \code{DNAStringSet} object with the genome sequence
#' @param type Either 'spliced' or 'unspliced'
#'
#' @return A \code{DNAStringSet} object with intronic sequences
#'
extractTxSeqs <- function(gtf, genome, type = "spliced") {
  ## Construct TxDb from gtf file. 
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf, format = "gtf")
  
  ## Group exons by transcript. When using exonsBy with by = "tx", 
  ## the returned exons are ordered by ascending rank for each transcript, 
  ## that is, by their position in the transcript. 
  grl <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
  
  ## Extract transcript sequences.
  ## Here, it's important that for each transcript, the exons must be ordered 
  ## by ascending rank, that is, by ascending position in the transcript.
  if (type == "spliced") {
    txout <- GenomicFeatures::extractTranscriptSeqs(x = genome, transcripts = grl)
  } else if (type == "unspliced") {
    txout <- GenomicFeatures::extractTranscriptSeqs(x = genome, transcripts = range(grl))
  } else {
    stop("Unknown 'type' argument")
  }
  
  return(txout)
}

#' Extract intron sequences
#'
#' @param gtf The path to a gtf file
#' @param genome A \code{DNAStringSet} object with the genome sequence
#' @param type Either 'collapse' or 'separate'
#' @param flanklength The length of the exonic flanking sequence
#' @param joinOverlappingIntrons Whether overlapping intron sequences (after adding 
#'   the flanking sequence) should be joined into a single intron
#'
#' @return A \code{DNAStringSet} object with intronic sequences
#' 
extractIntronSeqs <- function(gtf, genome, type = "collapse", flanklength = 90,
                              joinOverlappingIntrons = FALSE) {
  ## Construct TxDb from gtf file
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf, format = "gtf")
  
  if (type == "separate") {
    ## Group exons by transcript
    grl <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
  } else if (type == "collapse") {
    ## Group exons by gene
    grl <- GenomicFeatures::exonsBy(txdb, by = "gene")
    
    ## Collapse the exons of each gene
    grl <- GenomicRanges::reduce(grl)
  } else {
    stop("Unknown 'type' argument")
  }
  
  ## Get introns as the set difference between the range and the exons,
  ## for each transcript/gene
  ## Here, the order of the introns doesn't really matter, since they
  ## will be considered separately (not joined together into a transcript)
  grl <- BiocGenerics::setdiff(range(grl), grl)
  
  ## Add flanking region
  grl <- grl + flanklength
  
  if (joinOverlappingIntrons) {
    ## If two (introns + flanklength) overlap, join them
    grl <- GenomicRanges::reduce(grl)
  }
  
  gr <- unlist(grl)
  
  ## Add -I{X} to names
  names(gr) <- gsub("\\-I\\.", "-I", make.unique(paste0(names(gr), "-I")))
  
  ## Get sequence
  gs <- BSgenome::getSeq(x = genome, names = gr)
  
  ## Manually set names of extracted sequences
  stopifnot(all(width(gs) == width(gr)))
  names(gs) <- names(gr)
  
  return(gs)
}