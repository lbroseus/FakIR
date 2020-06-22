# fakIR: First-aid kit for studying IR events using RNA-sequencing data
# Copyright (C) 2019-2020   <lucile.broseus@igh.cnrs.fr>
#########################################################
#
# Extract reads overlapping at least one intron interval
#
#########################################################

#_______________________________________________________#
#' Extract reads overlapping at least one intron interval
#'
#' @description This function extracts (long) reads overlapping at least one intron interval.
#' Intron interval should be either specified as a \code{GRanges} object in \code{intronGR}, or will
#' be defined internally if a gtf file is provided in \code{gtf}.
#'
#' NB: it takes into account possible gaps and splicing in read alignments.
#'
#' @import GenomicRanges
#' @import IRanges
#' @import GenomicFeatures
#'
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom Rsamtools ScanBamParam
#' @importFrom S4Vectors Hits queryHits subjectHits
#' @importFrom rtracklayer export
#' @importFrom GenomicAlignments readGappedReads
#'
#' @param alignments A path to long read alignments in a bam file or a \code{GappedReads} object with long read alignments.
#' @param intronGR A \code{GRanges} object with intronic interval.
#' @param gtf Path to a reference transcriptome in a gtf file from which to extract intron intervals.
#' @param saveFile Bam file where long reads overlapping an intron will be saved.
#' @param keepSecondaryAlignments Whether secondary alignments should be considered.
#'                                Only considered when \code{alignments} is a file. (Default: FALSE)
#' @param overlapType Type of intron overlap to consider.
#'                    Either "within" to select only read that contain at least one intron to its full-extent.
#'                    Or "any" to consider partial intron overlaps as well. (Default: "any")
#' @param intronMinoverlap Minimum overlap (in bp) between a read alignment and an intron. (Default: 5)
#' @param verbose Whether to display some figures along the filtering process. (Default: FALSE)
#'
#' @return Either returns a \code{GappedReads} object if \code{saveFile=NULL} or outputs a bam file
#' containing input long reads overlapping at least one intron into \code{saveFile} .
#'
#' @seealso ['createIntronGR']
#'
#' @export
#--------------------------------------------------------#
#

extractIRreads <- function(alignments,
                           intronGR = NULL, gtf = NULL,
                           saveFile = NULL,
                           intronMinoverlap = 5,
                           overlapType = "any",
                           keepSecondaryAlignments = FALSE,
                           verbose = FALSE){

  if( missing(gtf) & missing(intronGR) ) stop( "Either gtf or intronGR must be specified..." )

  if( !missing(gtf) & !missing(intronGR) ) stop( "gtf or intronGR should not be set simultaneously..." )

  if( overlapType %in% c("within", 'any') ) stop( "overlapType value is invalid. It must be 'within'' or 'any'." )

  ## If needed call createIntronGR() for intronic intervals:
  if( !is.null(gtf) ){

    if( verbose ) cat("[IRFindeR2] IR-extraction:", "defining introns from gtf:\n", gtf, "\n")
    intronGR <- createIntronGR(gtf = gtf)

  }

  ## If needed, import splice alignments:
  if( !is(object = alignments, class2 = "GappedReads") ){

    param <- ScanBamParam(scanBamFlag(isUnmappedQuery = F, isSecondaryAlignment = keepSecondaryAlignments))
    alignments <- readGappedReads(file = alignments, param = param, use.names = TRUE)

  }
  seqlevelsStyle(alignments) <- seqlevelsStyle(intronGR) <- "Ensembl"

  if( verbose ) cat("[IRFindeR2] IR-extraction:", "analysing", length(alignments), "entries \n")

  ## Find overlaps between introns and read alignments:
  hits <- findOverlaps(query = intronGR, subject = alignments, minoverlap = intronMinoverlap, type = "within")

  ## Subset alignments with intron overlap:
  alignments <- alignments[ unique(subjectHits(hits)) ]
  if( verbose ) cat( "[IRFindeR2] IR-extraction:", length(alignments), " alignments kept. \n")


  if( !is.null( saveFile ) ){
    export(object = alignments, con = saveFile, format = "bam")
    if( verbose ) cat( "[IRFindeR2] IR-extraction:", length(alignments), " alignments saved in: \n", saveFile, "\n")
  }
  else return( alignments )

}
