# fakIR: First-aid kit for studying IR events using RNA-sequencing data
# Copyright (C) 2019-2020   <lucile.broseus@igh.cnrs.fr>
#########################################################
# FUNCTIONS:
#
# filterIRtranscripts()
#-------------------------------------------------------#
#
#########################################################

#_______________________________________________________#
#' Extract assembled Intron-retaining transcripts  using
#' a reference annotation
#'
#' @importFrom rtracklayer import
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @importFrom GenomicRanges GRanges reduce intersect findOverlaps
#' @importFrom GenomicFeatures intronsByTranscript makeTxDbFromGRanges
#' @importFrom S4Vectors Hits queryHits subjectHits
#'
#' @param assembly.gtf Path to a (custom) transcript assembly in a gtf file
#' @param reference.gtf Path to an extensive reference transcript assembly in a gtf file from which introns will be defined
#' @param includeAnnotatedIR \code{Boolean}. Whether to consider introns retained in a reference transcript (Default: TRUE).
#'
#' @return \code{GRanges} or gtf containing intron-retaining transcripts in the assembly
#'
#' @export
#--------------------------------------------------------#
#

filterIRtranscripts <- function(reference.gtf, assembly.gtf, includeAnnotatedIR = TRUE){

  reference <- rtracklayer::import( reference.gtf )

  # Remove IR transcripts from reference
  if( includeAnnotatedIR ){

    col <- intersect(grep(x = colnames(mcols(reference)), pattern = "transcript"),
                     grep(x = colnames(mcols(reference)), pattern = "type"))[1]

    colnames(mcols(reference))[col] <- "transcript_type"

    reference <- reference[ -which(reference$transcript_type == "retained_intron") ]

  }

  ## List intronic parts
  introns <- suppressWarnings( intronsByTranscript(x = makeTxDbFromGRanges(reference)) )
  introns  <- reduce(introns)
  #rm( reference )

  assembly <- rtracklayer::import(assembly.gtf)
  exons <- assembly[ which(assembly$type == "exon") ]

  seqlevelsStyle(exons) <- seqlevelsStyle(introns) <- "Ensembl"

  hits <- findOverlaps(query = introns, subject = exons , type = "within")

  tx_ids <- exons$transcript_id[ subjectHits(hits) ]
  gene_ids <- exons$gene_id[ subjectHits(hits) ]

  assembly <- assembly[ which( (assembly$type == "gene" & assembly$gene_id %in% gene_ids) |
                                 (assembly$transcript_id %in% tx_ids) )  ]

  return( assembly )

}
