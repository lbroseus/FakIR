# fakIR: First-aid kit for studying IR events using RNA-sequencing data
# Copyright (C) 2019-2020   <lucile.broseus@igh.cnrs.fr>


#_______________________________________________________#
#' findCompatibleTranscripts()
#' Find transcripts compatible with detected intron retention events
#'
#' @description This function finds transcripts compatible with detected intron retention events
#'
#' @import magrittr
#' @import GenomicRanges
#' @import GenomicFeatures
#'
#' @importFrom utils globalVariables
#' @importFrom rtracklayer import export
#' @importFrom dplyr arrange mutate distinct
#' @importFrom GenomeInfoDb seqlevelsStyle
#'
#' @param ir_events Object that can be coerced to a \code{GRanges} object, containing intron genomic intervals
#' @param gtf Path to the reference transcriptome where to search for transcripts with intron retention events
#' @param transcript_biotypes specifies biotype(s) of annotated transcripts to consider (Default: NULL -ie: all biotypes)
#' @param outfile (Defaut: NULL)
#'
#' @return If a valid path is specified in \code{outfile}, saves hypothetical intron-retaining transcripts in gtf format,
#'         otherwise returns results in a \code{Ranges} object.
#'
#' @export
#--------------------------------------------------------#

findCompatibleTranscripts <- function(ir_events, gtf, transcript_biotypes = NULL, outfile = NULL){

  #globalVariables(c("transcript_biotype", "transcript_id"))

  ir_events <- GRanges( ir_events )
  ### Adjust intron coordinates
  start(ir_events) <- start(ir_events)+1

  tx <- import( gtf )

  ## Selection of Annotated Transcripts according to biotype
  if( ! is.null(transcript_biotypes) ){

    biotype_col <- intersect(grep(x = colnames(mcols(tx)), pattern = "transcript"),
                             grep(x = colnames(mcols(tx)), pattern = "type"))

    colnames(mcols(tx))[ biotype_col ] <- "transcript_biotype"

    tx <- tx[ which( tx$type == "gene" | tx$transcript_biotype %in% transcript_biotypes ) ]

    tx <- tx[ which( tx$type == "gene" | tx$transcript_id %in% tx$transcript_id[tx$type == "transcript"] ) ]

  }

  ## Determine introns contained in each transcript
  intron_by_tx <- suppressWarnings( intronsByTranscript(x = makeTxDbFromGRanges(tx), use.names = T) )
  intron_by_tx <- unlist(intron_by_tx)
  intron_by_tx$transcript_id <- names(intron_by_tx)
  names(intron_by_tx) <- NULL

  ## Find detected retained introns matching transcripts' intron(s)
  seqlevelsStyle(intron_by_tx) <- seqlevelsStyle(ir_events) <- "Ensembl"

  matches <- findOverlaps(query = ir_events, subject = intron_by_tx, type = "equal")

  matches <- cbind.data.frame(ir_events[queryHits(matches)],
                              transcript_id = intron_by_tx$transcript_id[ subjectHits(matches)])

  matches <- matches %>% mutate(intron_id = paste(seqnames, ":", start, "-", end, ":", strand, sep = ""))

  ### Assemble Retained Intron and compatible exonic structures
  tx_structure <- suppressWarnings( exonsBy(x = makeTxDbFromGRanges(tx), use.names = T, by = "tx") )

  ### List exonic structures parallel to detected introns
  matches$transcript_id <-  as.vector(matches$transcript_id)
  index <- tapply(X = matches$transcript_id,
                  INDEX = 1:nrow(matches),
                  FUN = function(y) grep(pattern = y, x = names(tx_structure)))
  index <- as.vector( index )
  tx_structure <- tx_structure[ index ]

  ## Get gene rows in gtf formats
  genes <- tx[ which(tx$type == "gene" & tx$gene_id %in% tx$gene_id[tx$transcript_id %in% names(tx_structure)]) ]

  ## Get transcripts rows in gtf formats (possibly duplicated)
  index <- tapply(X = names(tx_structure),
                  INDEX = 1:length(tx_structure),
                  FUN = function(y) grep(pattern = y, x = tx$transcript_id[ which(tx$type == "transcript") ]))
  index <- as.vector( index )
  transcripts <-  tx[ which(tx$type == "transcript")][ index ]

  ## Build hypothetical IR-transcripts' structures
  tx_structure <- punion(tx_structure, GRanges(matches))
  names(tx_structure) <- paste(names(tx_structure), "+(", matches$intron_id, ")", sep = "")

  ## Adapt annotation ids
  transcripts$transcript_id <- names(tx_structure)

  tx_structure <- unlist( tx_structure )
  tx_structure$transcript_id <- names( tx_structure )
  tx_structure <- tx_structure[ order(tx_structure$transcript_id) ]

  ## Fromat exons
  tx_structure$type <- "exon"
  names(tx_structure) <- NULL

  geneAndTx <- data.frame(gene_id = transcripts$gene_id, transcript_id = transcripts$transcript_id)
  geneAndTx <- geneAndTx %>% distinct(transcript_id, .keep_all = T)

  exonsGenes <- data.frame(transcript_id = tx_structure$transcript_id)
  exonsGenes <- merge(exonsGenes, geneAndTx, by = "transcript_id", all.x = T)
  exonsGenes <- arrange(exonsGenes, transcript_id)

  tx_structure$gene_id <- exonsGenes$gene_id

  tx <- c(genes, transcripts, tx_structure)
  tx <- c(transcripts, tx_structure)
  tx <- sort( tx )

  if( !is.null(outfile) ) export(tx, con = outfile)
  else return( tx )

}
