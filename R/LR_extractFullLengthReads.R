# fakIR: First-aid kit for studying IR events using RNA-sequencing data
# Copyright (C) 2019-2020   <lucile.broseus@igh.cnrs.fr>
#########################################################
#
# Extract full-length long reads using
# annotated transcripts or/and CAGE data
#
#########################################################

#_______________________________________________________#
#' Extract full-length long reads
#'
#' @description This function extracts long reads whose alignment start position falls nearby
#' a known Transcription Start Sites (TSS).
#' For so doing, it makes use of TSS extracted from a given reference transcriptome (annotated TSS)
#' and/or of TSS inferrred from CAGE data.
#' At least one gtf file (in \code{gtf} or \code{cageFile}) must be specified.
#'
#' At the moment, it is intended for direct RNA Nanopore data,
#' for which molecules are read from the 3' to the 5' end.
#'
#' NB: it takes into account possible gaps and splicing in read alignments.
#'
#' @details
#' Splice-aware genomics alignments are imported into R as a \code{GappedReads} object.
#' Long reads' alignment start position (most 3' end of the RNA molecule) is identified.
#' If a TSS from input \code{gtf} or \code{cageFile} files falls within a \code{distanceToTSS}-wide window
#' around the read start position, then the read is considered to be full-length.
#'
#' @import GenomicRanges
#' @import IRanges
#' @import GenomicFeatures
#' @import GenomeInfoDb
#' @import Rsamtools
#'
#' @importFrom S4Vectors Hits queryHits subjectHits
#' @importFrom rtracklayer import export
#'
#' @param alignments Path to long read alignments in a bam file or a \code{GappedReads} object with long read alignments
#' @param gtf Path to a reference transcriptome (gtf/gff file) from which to extract annotated TSS
#' @param cageFile Path to Transcription Start Sites inferred (eg: from CAGE data)
#' @param saveFile Name of the bam file that will be output
#' @param distanceToTSS Maximum distance to a TSS (in bp). (Default: 25)
#' @param keepSecondaryAlignment Whether secondary alignments should be considered. (Default: FALSE)
#' @param verbose Whether some information about read filtering should be displayed.(Default: FALSE)
#'
#' @return Either returns a \code{GappedReads} object if \code{saveFile=NULL} or outputs a bam file
#' containing read alignments of long reads that may be full-length .
#'
#' @examples
#' \dontrun{
#' # Long read splice-aware alignments on a reference genome (eg: output by Minimap2):
#' bamFile <- system.file("extdata", "MCF10A_chr2.dont.bam", package = "IRFindeR2")
#'
#' #Reference transcriptome in a gtf file:
#' gtf <- system.file("extdata", "Homo_sapiens.GRCh38.97.chr2.gtf", package = "IRFindeR2")
#'
#' # TSS called from CAGE data carried on a comparable cell line (here MCF7), provided by ENCODE:
#  cageFile <-  system.file("extdata", "GSM849364_hg38_wgEncodeRikenCageMcf7CellPapTssGencV7.chr2", package = "IRFindeR2")
#'
#' # Select full-length reads alignments
#' # using both annotated transcripts and called TSS (eg: from ENCODE):
#' full_length_reads <- extractFullLengthReads(bamFile, gtf, cageFile)
#'
#' # Select full-length reads alignments using only annotated transcripts:
#' full_length_reads <- extractFullLengthReads(bamFile, gtf)
#'
#' # Select full-length reads alignments using only called TSS (eg: from ENCODE or in-house data):
#' full_length_reads <- extractFullLengthReads(bamFile, cageFile)
#' }
#'
#' @export
#--------------------------------------------------------#
#

extractFullLengthReads <- function(alignments,
                                   gtf = NULL,
                                   cageFile = NULL,
                                   saveFile = NULL,
                                   distanceToTSS = 25,
                                   keepSecondaryAlignment = FALSE,
                                   verbose = FALSE){

  # Some checks
  if( missing(gtf) & missing(cageFile) ) stop( "Either gtf or cageFile must be specified..." )
  if( distanceToTSS < 1 ) stop( "distanceToTSS must be > 0..." )

  if( !is(object = alignments, class2 = "GappedReads") ){

    # Import (splice aware) alignments
    what <- c("flag", "cigar")
    flag <- scanBamFlag(isUnmappedQuery = F, isSecondaryAlignment = keepSecondaryAlignment)
    param <- ScanBamParam(what=what, flag=flag)
    #----------------------------------#
    alignments <- readGappedReads(file = alignments, param = param, use.names = TRUE)

  }
  seqlevelsStyle( alignments ) <- "Ensembl"

  if( verbose ) cat("[IRFindeR2] FL-extraction:", "analysing", length(alignments), "entries \n")

  readStartWindows <- GRanges(seqnames = seqnames(alignments),
                                strand = strand(alignments),
                                ranges = IRanges(start = ifelse(strand(alignments)=="+",
                                                            start(alignments)-distanceToTSS, end(alignments)-distanceToTSS),
                                             end =  ifelse(strand(alignments)=="+",
                                                           start(alignments)+distanceToTSS, end(alignments)+distanceToTSS)))

  indices <- c()

  # Import tx data
  if( !is.null(gtf) ){

    promoters <- import(con = gtf)
    seqlevelsStyle( promoters ) <- "Ensembl"

    promoters <- promoters(promoters, upstream = 1, downstream = 1)

    hits <- findOverlaps(query = promoters, subject = readStartWindows)
    indices <- union(indices, subjectHits(hits))

  }

  # Import cage data
  if( !is.null(cageFile) ){

    cage <- rtracklayer::import(con = cageFile)
    seqlevelsStyle( cage ) <- "Ensembl"

    hits <- findOverlaps(query = cage, subject = readStartWindows)
    indices <- union(indices, subjectHits(hits))

  }

  alignments <- alignments[ unique( indices ) ]

  if( verbose ) cat( "[IRFindeR2] FL-extraction:", length(alignments), " alignments kept. \n")


  if( !is.null(saveFile) ){
    export(object = alignments, con = saveFile, format = "bam")
    if( verbose ) cat( "[IRFindeR2] FL-extraction:", length(alignments), " alignments saved in: \n", saveFile, "\n")
  }
  else return( alignments )

}
