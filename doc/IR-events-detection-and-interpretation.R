## ---- include = FALSE, out.width=10-------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----git_install, eval = FALSE, echo = TRUE-----------------------------------
#  
#  install_gitlab("lbroseus/IRFindeR2", host = "gitlab.igh.cnrs.fr")
#  

## ----packages, echo = T, eval = T, warning=F----------------------------------

suppressPackageStartupMessages( library(IRFindeR2) )
 

## ----other_packages, echo = T, eval = T, warning=F----------------------------

if( !("dplyr" %in% installed.packages()) ) install.packages( "dplyr" )

suppressPackageStartupMessages( library(dplyr) )

if( !("ggplot2" %in% installed.packages()) ) install.packages( "ggplot2" )

suppressPackageStartupMessages( library(ggplot2) )
 

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  
#  system.file("extdata", package = "IRFindeR2")
#  

## ---- eval = TRUE, echo = TRUE------------------------------------------------

gtf <- system.file("extdata", "Homo_sapiens.GRCh38.97.chr2.gtf", package = "IRFindeR2")
 

## ---- eval = TRUE, echo = TRUE------------------------------------------------

cageFile <- system.file("extdata", "GSM849364_hg38_wgEncodeRikenCageMcf7CellPapTssGencV7.chr2.gtf", 
                       package = "IRFindeR2")
 


## ---- echo = TRUE, eval = TRUE------------------------------------------------

irfinderFile <- system.file("extdata", "MCF10A_chr2.IRFinder-IR-dir.txt", package = "IRFindeR2")
 

## ---- eval = TRUE, echo = TRUE------------------------------------------------

bamFile <- system.file("extdata", "MCF10A_chr2.dont.bam", package = "IRFindeR2")
 

## ---- eval = TRUE, echo = TRUE------------------------------------------------

assembly.gtf <- system.file("extdata", "MCF10A_chr2.dont_assembly.gtf", package = "IRFindeR2")
 

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  
#  # Display available genomes and identify the one that fits your data:
#  BSgenome::available.genomes()
#  
#  # Specify the reference of your choice:
#  yourGenomeChoice <- "BSgenome.Hsapiens.UCSC.hg38"
#  
#  # Download its representation in R:
#  if (interactive()) {
#      if (!require("BiocManager"))
#          install.packages("BiocManager")
#      BiocManager::install( yourGenomeChoice )
#  }
#  
#  # Import into R:
#  BSgenome <- getBSgenome(genome = yourGenomeChoice, masked=FALSE)
#  

## ---- echo = TRUE, eval = TRUE------------------------------------------------

BSgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
 

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  
#  #Step 1: Build a genome index:
#  
#  bin/IRFinder -m BuildRef -r REF/Human-hg38-release97 \
#      -e REF/extra-input-files/RNA.SpikeIn.ERCC.fasta.gz \
#      -R REF/extra-input-files/Human_hg38_nonPolyA_ROI.bed \
#      ftp://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf.gz
#  

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  
#  #Step 2: align (paired-end) reads with STAR and analyse intron coverage:
#  
#  bin/IRFinder -r REF/Human-hg38-release97  -d irfinder READ_1.fastq READ_2.fastq
#  

## ---- echo = TRUE, eval = TRUE------------------------------------------------

irfinderFile <- system.file("extdata", "MCF10A_chr2.IRFinder-IR-dir.txt", package = "IRFindeR2")
 

## ---- echo = TRUE, eval = TRUE------------------------------------------------

irfinder_results <- read.table(irfinderFile, header = T, sep = "\t")

head(irfinder_results, 10)
 

## ----ir_detection, echo = TRUE, eval = TRUE-----------------------------------

# Main variables from which IR events are inferred:
minIRratio <- 0.1
minSpliceRead <- 3
minExonToIntronReads <- 1
minCoverage <- 0.1

# Import the file IRFinder-IR-dir.txt and extract filter IR events:
ir_events <- filterIRevents(irfinderFile = irfinderFile, 
                            minIRratio = minIRratio, 
                            minSpliceRead = minSpliceRead, 
                            minExonToIntronReads = minExonToIntronReads,
                            minCoverage = minCoverage)
 

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  
#  gtf <- "Your_favourite_organism_versionXXX.gtf"
#  

## ---- eval = TRUE, echo = TRUE------------------------------------------------

files <- c("sample1.WT.bam", "sample2.WT.bam", "sample3.WT.bam",
           "sample1.KO.bam", "sample2.KO.bam", "sample3.KO.bam")

conditions <- c("WT", "WT", "WT", 
                "KO", "KO", "KO") 

libType <- "single-end" # choose between "single-end" or "paired-end"

target <- data.frame(file = files, condition = conditions)

target
  

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  
#  idxd <- DEXSeq4Introns(gtf, target, libType, verbose = TRUE)
#  

## ---- eval = FALSE, echo = TRUE-----------------------------------------------
#  
#  table ( idxd$padj < 0.1 )
#  

## ----make_introns_db, eval = TRUE, echo = TRUE--------------------------------

intronGR <- createIntronGR(gtf = gtf)
 

## ---- eval = TRUE, echo = TRUE------------------------------------------------

intronGR
 

## ---- eval = TRUE, echo = TRUE------------------------------------------------

table( intronGR$annotated_intron )
 

## ---- eval = TRUE, echo = TRUE, fig.width = 10, fig.height = 5----------------

data.frame(intronGR) %>% mutate(intron_width = end - start + 1) %>% 
                         ggplot(aes(intron_width, fill = annotated_intron)) + 
                         geom_histogram(bins = 25, color = "orange") + scale_x_log10() +
                         theme_minimal() + scale_fill_manual(values = c("darkblue", "lightblue")) +
                         ggtitle("Independent intron width (Homo Sapiens, chromosome 2)")
 

## ----call_events, eval = TRUE, echo = TRUE------------------------------------

IR_events <- callIRevents(bamFile, intronGR)
 

## ----display_events, eval = TRUE, echo = TRUE---------------------------------

head(IR_events)
 

## ----full_ir_events, eval = TRUE, echo = TRUE---------------------------------

full_IR_events <- IR_events %>% 
                  filter(intron_fracoverlap == 1 & intron_startOverlap & intron_endOverlap)
 

## ----reliable_events, eval = TRUE, echo = TRUE--------------------------------

retained_introns <- IR_events %>% 
                    filter(intron_fracoverlap == 1 & intron_startOverlap & intron_endOverlap) %>%
                    group_by(intron_chr, intron_start, intron_end, intron_strand, gene_id) %>%
                    summarise(read_count = n()) %>%
                    data.frame()
 

## ----IR_rna, eval = TRUE, echo = TRUE-----------------------------------------

IR_reads <- full_IR_events %>% 
            group_by(read_name) %>% 
            summarise(nbfullRI = n()) %>% 
            data.frame()
 

## ---- eval = F, echo = F------------------------------------------------------
#  
#  #Multiple retention events
#  
#  #Long read depth analysis
#  

## ----stringtieA, eval = FALSE, echo = TRUE------------------------------------
#  
#  stringtie -L -G Homo_sapiens.GRCh38.97.gtf -o assembly.gtf long_reads.bam
#  

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  
#  # Long read splice-aware alignments on a reference genome (eg: output by Minimap2):
#  alignments <- system.file("extdata", "MCF10A_chr2.dont.bam", package = "IRFindeR2")
#  
#  # Reference transcriptome in a gtf file:
#  gtf <- system.file("extdata", "Homo_sapiens.GRCh38.97.chr2.gtf", package = "IRFindeR2")
#  
#  # TSS called from CAGE data carried on a comparable cell line (here MCF7), provided by ENCODE:
#  cageFile <-  system.file("extdata", "GSM849364_hg38_wgEncodeRikenCageMcf7CellPapTssGencV7.chr2.gtf",
#                           package = "IRFindeR2")
#  
#  # Name for the output bam file, that will contain only full-length reads alignments:
#  saveFile <- "fl_long_reads.bam"  # Or replace with a convenient path
#  

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  
#  extractFullLengthReads(alignments, gtf, cageFile, saveFile, verbose = TRUE)
#  

## ---- eval = TRUE, echo = TRUE------------------------------------------------

gtf <- system.file("extdata", "Homo_sapiens.GRCh38.97.chr2.gtf", package = "IRFindeR2")

bamFile <- system.file("extdata", "MCF10A_chr2.dont.bam", package = "IRFindeR2")

saveFile <- "ir_long_reads.bam"
# Or replace with a convenient path


## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  
#  extractIRreads(alignments, gtf = gtf, saveFile = saveFile, intronMinoverlap = 5, keepSecondaryAlignments = FALSE)
#  

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  
#  # Define intronic intervals from a reference transcriptome
#  # Or feed intronGR with custom intron interval in a GRanges
#  intronGR <- createIntronGR( gtf )
#  
#  extractIRreads(alignments, intronGR = intronGR, saveFile, intronMinoverlap = 5, keepSecondaryAlignments = FALSE)
#  

## ----stringtieIR, eval = FALSE, echo = TRUE-----------------------------------
#  
#  stringtie -L -G Homo_sapiens.GRCh38.97.gtf -f 0 -o assembly.gtf ir_long_reads.bam
#  

## ----genome, echo = TRUE, eval = TRUE-----------------------------------------

BSgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
 

## ----stopCodons, echo = TRUE, eval = TRUE-------------------------------------

ir_events <- countStopCodons(ir_events, BSgenome, stopCodons = c("TAA", "TAG", "TGA"))
 

## ---- echo = TRUE, eval = TRUE------------------------------------------------

ir_events %>% dplyr::select(seqnames, start, end, strand , frame1, frame2, frame3) %>% head(15)
 

## ----NMDplus, echo = TRUE, eval = TRUE----------------------------------------

ir_events %>% filter(frame1*frame2*frame3 > 0) %>% nrow()

ir_events %>% filter(frame1*frame2*frame3 > 0) %>% 
              dplyr::select(seqnames, start, end, start, gene_name, frame1, frame2, frame3) %>%
              head(15)
 

## ----NMDminus, echo = TRUE, eval = FALSE--------------------------------------
#  
#  ir_events %>% filter(frame1+frame2+frame3 == 0)
#  

## ----NMD0, echo = TRUE, eval = TRUE-------------------------------------------

ir_events %>% filter(frame1+frame2+frame3 != 0 & frame1*frame2*frame3 == 0) %>%
              dplyr::select(seqnames, start, end, start, gene_id, frame1, frame2, frame3)
 

## ----ref, echo = TRUE, eval = TRUE--------------------------------------------

gtf <- system.file("extdata", "Homo_sapiens.GRCh38.97.chr2.gtf", package = "IRFindeR2")
 

## ----assembly, echo = TRUE, eval = TRUE---------------------------------------

# Where to save IR-transcripts assembly:
outfile <- paste(system.file("extdata",  package = "IRFindeR2"), "MCF10A.hypothetical_ir_transcripts.gtf", sep = "/")

findCompatibleTranscripts(ir_events, gtf, transcript_biotypes = "protein_coding", outfile = outfile)
 

## ----isa_install, eval = FALSE, echo = TRUE-----------------------------------
#  
#  install.packages("BiocManager")
#  BiocManager::install()
#  
#  BiocManager::install("IsoformSwitchAnalyzeR")
#  

## ----isa, eval = T, echo = T, warning = F-------------------------------------

suppressPackageStartupMessages( require(IsoformSwitchAnalyzeR) )
 

## ---- eval = T, echo = T------------------------------------------------------

assembly.gtf <- paste(system.file("extdata",  package = "IRFindeR2"), "MCF10A_chr2.dont_assembly.gtf", sep = "/")
  

## ----importGTF, eval = T, echo = T--------------------------------------------

mySwitchList <- IsoformSwitchAnalyzeR::importGTF(pathToGTF = assembly.gtf)
 

## ----genome2, eval = T, echo = T----------------------------------------------

BSgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38


## ----isa_vignette, eval = T, echo = T, warning = F----------------------------

vignette( "IsoformSwitchAnalyzeR" )
 

## ----ORFs, eval = T, echo = T-------------------------------------------------

orfMethod <- "mostUpstream"

mySwitchList <- analyzeORF(switchAnalyzeRlist = mySwitchList, 
                           genomeObject = BSgenome,
                           orfMethod = orfMethod, 
                           quiet = TRUE)
 

## ---- eval = F, echo = T------------------------------------------------------
#  
#  mySwitchList$orfAnalysis
#  

## ---- eval = T, echo = T------------------------------------------------------

mySwitchList$orfAnalysis %>% dplyr::select(isoform_id, 
                                           orfStartGenomic, orfEndGenomic, 
                                           stopDistanceToLastJunction, PTC) %>%
                             head(15)
 

## ---- eval = T, echo = T------------------------------------------------------

table( mySwitchList$orfAnalysis$PTC )


