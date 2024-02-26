## #############################################################################
## Date:        March 023
## Author:      Allison M. Burns
## Filename:    2_PeakCalling_csaw.R
## Project:     Cas9 ChiP
## Description: Read in bam files and check fragment sizes, define enriched
##              windows and get counts, normalize using both TMM and Loess
##              methods for comparison.
## #############################################################################

library(csaw)
library(edgeR)
library(tidyverse)

## Load bam file information
bam.fls <- list.files("./data/alignment/2_mark_dups", pattern = ".bam$", full.names=TRUE)
bam.names <- substr(gsub("markdups.bam", "", basename(bam.fls)), 7, 26)
bam.names <- gsub("\\_$","", bam.names)
names(bam.fls) <- bam.names
bam.desc <- substr(bam.names, 1, nchar(bam.names)-5)

################################################################################
## 1. Reading bam files
################################################################################
## Set counting parameters
param <- readParam(pe = "both", ## account for paired end data
                   max.frag = 500, ## run get PEsizes to help decide this
                   dedup = TRUE, ## ignore marked duplicates
                   ## restrict = "15",
                   minq=20 ## minimum mapping quality
                   )

## Get paired-end fragment size distribution
frag.len <- lapply(bam.fls, function(bam) { getPESizes(bam, param = param) })
names(frag.len) <- bam.names

diagnostic <- do.call(rbind, lapply(frag.len, function(out) {
    c(out$diagnostics, too.large=sum(out$sizes > 1000))
}))

frag.size <- do.call(rbind,lapply(seq(1:length(frag.len)), function(i) {
    nom <- names(frag.len[i])
    data.frame(size = frag.len[[i]]$sizes,
               nom)
}))

## Plot 1/5th of the data (takes way too long otherwise
tt <- frag.size[sample(seq(1:nrow(frag.size)), nrow(frag.size)/5),]

pdf("./figures/alignment/bam_file_fragsizes.pdf", height = 5, width = 10)
ggplot(tt, aes(x = size, fill = nom)) +
    geom_histogram(binwidth=10) +
    facet_wrap(~nom) +
    theme_bw() +
    xlim(0,1000)
dev.off()

################################################################################
## 2. Get csaw denovo windows and read counts
################################################################################
## Choose small window size to account for small cas9 binding region
win.width <- 10
data <- windowCounts(bam.fls, width=win.width, param=param)

## filter lowly expressed windows by local enrichment (local background: 2kb)
neighbor <- suppressWarnings(resize(rowRanges(data), width=2000, fix="center"))
wider <- regionCounts(bam.fls, regions=neighbor, param=param) # count reads in neighborhoods
filter.stat <- filterWindowsLocal(data, wider)
csaw.local.filt <- data[filter.stat$filter > log2(3),] # 3x enrichment over 2kb neighborhood abundance
rtracklayer::export.bed(rowRanges(csaw.local.filt), con = "./data/local_filt.bed")

################################################################################
## 3. Normalization
################################################################################
## count BAM background bins (for TMM normalization)
binned <- windowCounts(bam.fls, bin=TRUE, width=10000, param=param)

## TMM normalization based on binned counts
csaw.local.tmm <- normFactors(binned, se.out=csaw.local.filt)

## csaw loess-normalization
csaw.local.loess <- normOffsets(csaw.local.filt, se.out=TRUE) 

################################################################################
## Save analysis
################################################################################
saveRDS(binned, "./data/deAnalysis/binned.rds")
saveRDS(csaw.local.tmm, "./data/deAnalysis/csaw_local_tmm.rds")
saveRDS(csaw.local.loess, "./data/deAnalysis/csaw_local_loess.rds")
