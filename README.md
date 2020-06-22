# FakIR: a First-aid kit for the computational prediction of Intron Retention events
===========================================================

FakIR is an R package for Intron Retention (IR) detection, quantification and exploration using Second and Third Generation Sequencing data.

A pdf version of the tutorial is available here: [Tutorial](https://github.com/lbroseus/FakIR/blob/master/IR-events-detection-and-interpretation.pdf).


## Package installation 

Please note that FakIR requires rlang version >= 0.4.5.

Also, FakIR makes use of functionalities from several packages from [Bioconductor](https://bioconductor.org) you might need to install first.
This can be done by copy-pasting the following code in a R session:

```

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install()


BiocManager::install(c("IsoformSwitchAnalyzeR", "DEXSeq", "BSgenome.Hsapiens.UCSC.hg38"))
  
```

Also, so as to install FakIR directly in a R session from GitHub, you will need to have *devtools* installed:

```

install.packages("devtools")
   
```

Then, install FakIR using:

```

devtools::install_github("lbroseus/FakIR", build_vignettes = T)
   
```

Alternatively, source files can be downloaded from the following link, and installed manually:

```

git clone https://github.com/lbroseus/FakIR.git
   
```

## Detecting Intron Retention events using RNA-seq data 

Typical workflows for detecting and analysing intron retention using RNA-seq data are detailed in the vignette and can be accessed in a R session using:

```

require(FakIR)

vignette("IR-events-detection-and-interpretation")
  
```

Especially, it illustrates how to perform the following tasks:

1. IR events detection from short read RNA-seq data
2. Differential analysis from short read RNA-seq data (using DEXSeq)
3. IR events detection from long read RNA-seq data
4. IR-transcripts reconstruction from long read RNA-seq data
5. Prediction of IR-NMD-targets  


## References

1. _Broseus and Ritchie **Challenges in detecting and quantifying intron retention from NGS data.** CSBJ. (2020)_

