# fakIR: First-aid kit for studying IR events using RNA-sequencing data
# Copyright (C) 2019-2020   <lucile.broseus@igh.cnrs.fr>
#########################################################
# FUNCTIONS:
#
# DEXSeq4IR()
#-------------------------------------------------------#
#
#########################################################


#_______________________________________________________#
#' Differential analysis of intronic counts using DEXSeq
#'
#'
#' @description This function is a wrapper adapting the DEXSeq (exon) method for differential
#' analysis to retained introns.
#'
#' @import magrittr
#' @import DEXSeq
#' @import GenomicFeatures
#' @import GenomicAlignments
#' @import Rsamtools
#'
#' @importFrom rtracklayer import
#' @importFrom S4Vectors Rle
#' @importFrom SummarizedExperiment colData
#' @importFrom DESeq2 plotDispEsts
#' @importFrom grDevices dev.off pdf
#'
#' @param gtf Path to a gtf file from which to extract (independent) intron intervals.
#' @param target A formatted \code{data.frame} summarising the sequencing experiment (with columns file and condition).
#' @param libType Library types of samples from the experiment. Either "single-end" or "paired-end".
#' @param keepAnnotatedRI Whether to analyse introns annotated as retained (Default: TRUE).
#'
#' @param plotDisp Whether to plot the graph of dispersion on screen (Default: TRUE).
#' @param saveDisp A path where to save the plot of dispersions (Default: NULL, no plot will be saved)
#' @param verbose Whether to output messages throughout the process (Default: FALSE).
#' @param noWarnings Whether to silence possible warnings (Default: TRUE).
#'
#' @return A \code{DEXSeqResults} object as output by function ['DEXSeqResults'] from the DEXSeq package.
#'
#' @seealso ['makeTarget']
#'
#' @export
#'
#--------------------------------------------------------#
# TESTED

DEXSeq4Introns <- function(gtf, target, libType,
                           keepAnnotatedRI = TRUE,
                           plotDisp = TRUE, saveDisp = NULL,
                           verbose = FALSE, noWarnings = TRUE){

  if( !("file" %in% colnames(target)) ) stop("target must be a data.frame a column named \"file\". \n")
  if( !("condition" %in% colnames(target)) ) stop("target must be a data.frame a column named \"condition\". \n")

  if( !(libType %in% c("single-end", "paired-end")) ) stop("Please indicate a valid library type. \n")


  if( verbose ) cat("Extracting intronic parts...")
  gtf <- import(gtf)

  if( keepAnnotatedRI ){
    coln <- intersect(grep("transcript", x = colnames(mcols(gtf))), grep("type", x = colnames(mcols(gtf))))
    if( length(coln) == 1 ){
      colnames(mcols(gtf))[coln] <- "transcript_type"
      coln <- which(gtf$transcript_type == "retained_intron")
      if( length(coln) > 0 ) gtf <- gtf[ -coln ]
    }
  }

  intronicParts <- suppressWarnings(
                          intronicParts(txdb = makeTxDbFromGRanges(gtf),
                                        linked.to.single.gene.only = TRUE)
                                   )

  rm( gtf )

  if( verbose ) cat("OK. \n")

  bam_ls <- BamFileList(as.character(target$file), index = character(), yieldSize = 100000, obeyQname = TRUE )

  #Some trick
  orderByGeneName <- order(unlist(intronicParts$gene_id, use.names=FALSE))
  intronic_rle <- runLength(Rle(unlist(intronicParts$gene_id[ orderByGeneName ],use.names=FALSE)))
  intronicParts <- intronicParts[ orderByGeneName ]
  intronic_part <- unlist(lapply(intronic_rle, seq_len), use.names=FALSE)
  intronicParts$exonic_part <- intronic_part

  if( verbose ) cat("Counting intronic reads...may take some time...")

  if( noWarnings) SE <- suppressWarnings( summarizeOverlaps( intronicParts, bam_ls, mode = "Union",
                           singleEnd = (libType == "single-end"),
                           ignore.strand = FALSE, inter.feature = FALSE, fragments = FALSE) )
  else SE <- summarizeOverlaps( intronicParts, bam_ls, mode = "Union",
                                                  singleEnd = (libType == "single-end"),
                                                  ignore.strand = FALSE, inter.feature = FALSE, fragments = FALSE)
  if( verbose ) cat("OK. \n")

  colData(SE)$condition <- target$condition

  if( verbose ) cat("Now, proceeding with DEXseq...")
  if( noWarnings ) dxd <- suppressWarnings( DEXSeqDataSetFromSE(SE, design= ~ sample + exon + condition:exon) )
  else dxd <- DEXSeqDataSetFromSE(SE, design= ~ sample + exon + condition:exon)

  if( noWarnings ) suppressWarnings( dxd <- estimateSizeFactors( dxd ) )
  else dxd <- estimateSizeFactors( dxd )

  ## Estimation and visualisarion of dispersion values

  dxd <- estimateDispersions( dxd )

  if( plotDisp ) plotDispEsts(dxd, main = "DEXSeq-Intron: dispersion estimates")

  if( !is.null(saveDisp) ){

    pdf( saveDisp )
    plotDispEsts(dxd, main = "DEXseq-IR: dispersion estimates")
    dev.off()

  }

  dxd <- testForDEU( dxd )

  if( noWarnings ) suppressWarnings( dxd <- estimateExonFoldChanges( dxd, fitExpToVar = "condition") )
  else dxd <- estimateExonFoldChanges( dxd, fitExpToVar = "condition")

  if( verbose ) cat("OK. \n")

  dxd <- DEXSeqResults( dxd )

  return( dxd )
}
