# fakIR: First-aid kit for studying IR events using RNA-sequencing data
# Copyright (C) 2019-2020   <lucile.broseus@igh.cnrs.fr>
#########################################################
# FUNCTIONS:
#
# filterIR()
#-------------------------------------------------------#
#
#########################################################

#_______________________________________________________#
#' Call Intron Retention (IR) events captured by NGS data using IRFinder results
#'
#'
#' @description This function is intended to call Intron Retention events from IRFinder analysis by applying
#' advised filters to a IRFinder-IR-dir.txt file.
#'
#' @import magrittr
#'
#' @importFrom dplyr filter select
#' @importFrom data.table fread
#' @importFrom stringr str_split
#' @importFrom utils globalVariables
#'
#' @param irfinderFile Path towards a IRFinder-IR-dir.txt file output by IRFinder
#' @param minIRratio Minimum retention level (\code{IRratio}) for an intron (value between 0 and 1)
#' @param minSpliceRead Minimum number of splice reads (\code{SpliceExact}) supporting the intron (Default value: 3)
#' @param minExonToIntronReads Minimum number of reads spanning the exon-to-intron junctions on both side of the intron (Default value: 1)
#' @param minCoverage Minimum fraction of the intron that should be covered by reads (Default value: 0.5)
#' @param onlyCleanIntrons Whether only introns without any known overlapping event (eg: exon) should be kept (Default: TRUE)
#'
#' @return A \code{data.frame} derived from \code{irfinderFile} with filtered IR events.
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#'  # Import formatted data.frame from IRfinder output, without any selection:
#'  filterIRevents(irfinderFile = irfinderFile,
#'                 minIRratio = 0, minSpliceRead = 0,
#'                 minExonToIntronReads = 0, minCoverage = 0)
#'  # Import formatted data.frame from IRfinder output, without usual intron selection:
#'  filterIRevents(irfinderFile = irfinderFile)
#'
#' }
#--------------------------------------------------------#
# TESTED

filterIRevents <- function(irfinderFile,
                           minIRratio,
                           minSpliceRead = 3,
                           minExonToIntronReads = 1,
                           minCoverage = 0.5,
                           onlyCleanIntrons = TRUE){

  #globalVariables(c("SpliceExact", "ExonToIntronReadsLeft", "ExonToIntronReadsRight", "IRratio",
  #                  "Name", "Warnings", "intron_status"))

  ir_events <- fread( irfinderFile, data.table = F)

  infos <- tapply(ir_events$Name,
                  1:nrow(ir_events),
                  FUN = function(x) as.vector( str_split(string = x, pattern = "/")[[1]]) )
  infos <- do.call(rbind.data.frame, infos)
  colnames(infos) <- c("gene_name", "gene_id", "intron_status")

  ir_events <- cbind.data.frame(ir_events, infos); rm(infos)
  ir_events <- ir_events %>% select(-Name, -Warnings, -Null)

  if( onlyCleanIntrons ) ir_events <- ir_events %>% filter(intron_status == "clean")

  ir_events <- ir_events %>%
               filter(SpliceExact >= minSpliceRead) %>%
               filter(ExonToIntronReadsLeft >= minExonToIntronReads & ExonToIntronReadsRight >= minExonToIntronReads) %>%
               filter(IRratio >= minIRratio)

  return( ir_events )

}
