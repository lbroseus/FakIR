# fakIR: First-aid kit for studying IR events using RNA-sequencing data
# Copyright (C) 2019-2020   <lucile.broseus@igh.cnrs.fr>
#########################################################
# FUNCTIONS:
#
# countStopCodons()
#-------------------------------------------------------#
#
#########################################################


#_______________________________________________________#
#' Count the number of codons in a sequence
#'
#' @description This function screens genomic interval (eg: intron) for stop codons
#'
#' @importFrom stringr str_length str_sub str_detect
#'
#' @param string A \cite{string} in which a set of codon sequences should be looked for
#' @param whichCodons A vector of strings specifying (stop) codons.
#'
#' @return A one-line \code{data.frame} giving the number of stop codons fround in each reading frame.
#'
#--------------------------------------------------------#

detectCodons <- function(string, whichCodons){

  nbs <- rep(0, 3)
  start_positions <- seq(from = 1, to = str_length(string)-2, by = 1)

  codons <- tapply(start_positions,
                   1:length(start_positions),
                   FUN = function(x) str_sub(string, start = x , end = x+2))

  for( f in 1:3 ){
    start_positions <- which(1:length(codons) %% f == 0)
    nbs[f] <- sum( tapply(whichCodons, 1:length(whichCodons),
                          FUN = function(x) sum(str_detect(string = codons[start_positions], pattern = x)) ) )
  }

  nbs <- data.frame( t(nbs) )
  colnames(nbs) <- paste("frame", 1:3, sep = "")

  return( nbs )

}

#_______________________________________________________#
#' Screen genomic intervals and count stop codons
#'
#' @description This function screens genomic interval (eg: intron) for stop codons
#'
#' @import GenomicRanges
#' @import GenomicFeatures
#'
#' @importFrom Biostrings getSeq
#' @importFrom GenomeInfoDb seqlevelsStyle
#'
#' @param granges An object (eg: a \code{data.frame}) than can be coerced to a \code{GRanges} object
#'                containing genomic intervals where stop codons will be searched for
#' @param BSgenome A \code{BSgenomeObject} specifying the reference genome sequence
#' @param stopCodons Vector of sequences defining DNA stop codons (Default: c("TAA", "TAG", "TGA"))
#'
#' @return A \code{data.frame} copy of \code{granges} to which three columns with the number of stop codons
#'         found in each reading frame are appended.
#'
#' @export
#'
#--------------------------------------------------------#

countStopCodons <- function(granges, BSgenome, stopCodons = c("TAA", "TAG", "TGA")){

  granges <- GRanges( granges )
  seqlevelsStyle(BSgenome) <- seqlevelsStyle(granges)[1]

  sequences <- getSeq(BSgenome, granges)
  sequences <- as.character( sequences )

  codon_detection <- tapply(sequences, 1:length(sequences),  function(x) detectCodons(x, whichCodons = stopCodons))
  codon_detection <- do.call(rbind.data.frame, codon_detection)

  granges <- cbind.data.frame(granges, codon_detection)

  return( granges )

}
