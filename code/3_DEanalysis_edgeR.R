## #############################################################################
## Date:        March 023
## Author:      Allison M. Burns
## Filename:    3_DEanalysis_edgeR.R
## Project:     Cas9 ChiP
## Description: Read in csaw defined peaks and run differential accessibility
##              analysis between IP and input groups, and merge nearby windows
##              to get DE peaks. Create volcano plot for DE peak regions and
##              assign peak regions to nearest promoter
## #############################################################################

library(csaw)
library(edgeR)
library(tidyverse)
library(GenomicFeatures)
library(ggrepel)
library(GenomicRanges)

## Load normalized window sets
comparisons <- list.files("./data/deAnalysis", full.names = TRUE, pattern = "rds")
names(comparisons) <- gsub(".rds","", basename(comparisons))
comps <- lapply(comparisons, function(x) { readRDS(x) })

binned <- comps[[1]]
comps <- comps[grep("^csaw",names(comps))]

## Load bam file information
bam.names <- colnames(comps[[1]])
bam.desc <- substr(bam.names, 1, nchar(bam.names)-5)

################################################################################
## DIFFERENTIAL ACCESSIBILITY ANALYSIS
################################################################################
final.merged.peaks <- lapply(seq(1:length(comps)), function(i) {
    nom <- names(comps[i])
    working.windows <- comps[[i]]
    
    ## setup design matrix
    y <- asDGEList(working.windows)
    colnames(y$counts) <- bam.names
    rownames(y$samples) <- bam.names
    y$samples$group <- bam.desc
    design <- model.matrix(~0+group, data=y$samples)
    colnames(design) <- unique(bam.desc)
    
    ## stabilize dispersion estimates with empirical bayes
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust=TRUE)
    
    ## Get replicate similarity
    pdf(paste("./figures/MDS_", nom, ".pdf", sep = ""), width =5, height = 5)
    adj.counts <- cpm(y, log=TRUE)
    plotMDS(adj.counts, main="top 5000 genes", col=c("blue", "blue", "red", "red"),
            labels=bam.names, top=5000)
    dev.off()
    
    ## testing for differentially-accessible windows
    arc.res <- glmQLFTest(fit,
                          contrast=makeContrasts(Arc_guide_flag-Arc_guide_input,
                                                 levels=design))
    nt.res <- glmQLFTest(fit,
                         contrast=makeContrasts(NT_guide_flag-NT_guide_input,
                                                levels=design))
    
    ## format DE results
    arc.windows <- working.windows
    rowData(arc.windows) <- arc.res$table
    nt.windows <- working.windows
    rowData(nt.windows) <- nt.res$table
    
    ## merge nearby windows if within 500 bp 
    arc.merged <- mergeWindows(rowRanges(arc.windows), tol=500L, max.width=5000L)
    nt.merged <- mergeWindows(rowRanges(nt.windows), tol=500L, max.width=5000L)
    
    ## use most significant window as stats representation for p-value and FDR for merged windows
    tab.best.arc <- getBestTest(arc.merged$id, arc.res$table)
    tab.best.nt <- getBestTest(nt.merged$id, nt.res$table)
    
    ## cat all relevant statistical data for final merged windows (no redundant columns)
    colnames(nt.res$table) <- paste(colnames(nt.res$table) ,"nt",sep = "_")
    
    final.merged.peaks <- GRanges(cbind(as.data.frame(arc.merged$region),
                                        arc.res$table[tab.best.arc$rep.test, -4],
                                        tab.best.arc[,-c(7:8)],
                                        nt.logFC = tab.best.nt$rep.logFC,
                                        nt.PValue = tab.best.nt$PValue,
                                        nt.FDR = tab.best.nt$FDR))
    
    ## sort by FDR (arc)
    final.merged.peaks <- final.merged.peaks[order(final.merged.peaks@elementMetadata$FDR), ]
    
    ## Save data
    saveRDS(final.merged.peaks,
            file = paste("./data/deAnalysis/DEpeaks_", nom, ".rds", sep = ""))
    writexl::write_xlsx(data.frame(final.merged.peaks),
                        path = paste("./data/deAnalysis/DEpeaks_", nom, ".xlsx", sep = ""))
    
    final.merged.peaks
})


################################################################################
## Add genes to excel files
################################################################################
## Load EdgeR peaks
de_peaks <- list.files("../data/deAnalysis", pattern = "DEpeaks", full.names = TRUE)
de_peaks <- de_peaks[grep(".rds", de_peaks)]
names(de_peaks) <- gsub("DEpeaks_|.rds","",basename(de_peaks))
de_peaks <- lapply(de_peaks, readRDS)
names(de_peaks) <- gsub("peak_counts", "macs2",names(de_peaks))
names(de_peaks) <- gsub("counts_local", "csaw",names(de_peaks))

## make txdb for mm10
## txdb <- makeTxDbFromEnsembl(organism = "Mus musculus", release = 102)
txdb <- loadDb("./data/annotations/mmusculus_ensembl_mm10_release102.txdb")

## Get promoters
proms <- promoters(txdb)

## Get common names
mm <- gprofiler2::gconvert(names(proms), organism= "mmusculus")
names(proms) <-  mm$name[match(names(proms), mm$input)]

## Find nearest gene
de_peak_labels <- lapply(de_peaks, function(x) {
    ## x <- de_peaks[[1]] ## for testing
    ## Find Closeset promoters
    x$closest_proms <- NA
    x$dist_closest_proms <- 0
    dist2prom <- distanceToNearest(x,proms)
    x[queryHits(dist2prom)]$closest_proms <- names(proms)[subjectHits(dist2prom)]
    x[queryHits(dist2prom)]$dist_closest_proms <- values(dist2prom)$distance
    x
})

saveRDS(de_peak_labels, file = "./data/deAnalysis/DEpeaks_labelled.rds")
writexl::write_xlsx(lapply(de_peak_labels, data.frame), ##data.frame(final.merged.peaks),
                    path = "./data/deAnalysis/DEpeaks_labelled.xlsx")

################################################################################
## Volcano Plots
################################################################################
peaks <- readRDS("./data/deAnalysis/DEpeaks_labelled.rds")

peaks.loess <- peaks[[2]]

peaks.arc <- data.frame(gn     = peaks.loess$closest_proms,
                        logFC  = peaks.loess$logFC,
                        PValue = peaks.loess$PValue,
                        FDR    = peaks.loess$FDR,
                        exp    = "Arc")
peaks.nt <- data.frame(gn     = peaks.loess$closest_proms,
                       logFC  = peaks.loess$nt.logFC,
                       PValue = peaks.loess$nt.PValue,
                       FDR    = peaks.loess$nt.FDR,
                       exp    = "NT")

peaks <- rbind(peaks.arc, peaks.nt)
peaks$sig <- case_when(peaks$logFC > 4 & peaks$FDR <= 0.05 ~ "up")

pdf("./figures/volcano_plot_local_loess.pdf", height = 5, width = 10)
ggplot(peaks, aes(x = logFC, y = -log10(PValue), color = sig, label = gn)) + 
    geom_point() +
    facet_wrap(~exp) +
    theme_bw() +
    geom_label_repel(data = subset(peaks, sig == "up"),
                     size = 3,
                     segment.size = 0.2)
dev.off()

################################################################################
## Define peak locations
################################################################################
## Load data
peaks <- readRDS("./data/deAnalysis/DEpeaks_labelled.rds")
peaks.loess <- peaks[[grep("loess", names(peaks))]]
## remove peaks that don't come from canonical chromosomes
peaks.loess <- peaks.loess[-grep("\\.",as.character(seqnames(peaks.loess)))] 

gns <- genes(txdb)
exons <- exonsBy(txdb, by = "gene")
proms <- promoters(txdb)

## Define Regions
peaks.loess$reg.type <- "intergenic"
peaks.loess$reg.type[queryHits(findOverlaps(peaks.loess, gns))] <-  "intronic"
peaks.loess$reg.type[queryHits(findOverlaps(peaks.loess, exons))] <-  "exonic"
peaks.loess$reg.type[queryHits(findOverlaps(peaks.loess, proms))] <-  "promoter"

## Set plotting color
cols <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")
names(cols) <- unique(peaks.loess$reg.type)

## Plot pie chart for ALL PEAKS
## Format dataframe for plotting
all.type <- data.frame(table(peaks.loess$reg.type))
all.type <- all.type |> 
    arrange(desc(Freq)) |>
    mutate(perc = round(100 * (Freq/sum(Freq)), 2)) 

all.peaks <- ggplot(all.type, aes(x="", y=perc, fill=Var1)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0) +
    theme_void() +
    geom_text(aes(label = paste(perc, "%", sep = "")),
              position = position_stack(vjust = 0.5),
              color = "white", size=2) +
    scale_fill_manual(values=cols) +
    labs(fill = "genomic region")

## Arc peaks
arc.loess <- peaks.loess[peaks.loess$FDR <= 0.05 & peaks.loess$logFC >= 4]
## Format dataframe for plotting
arc.type <- data.frame(table(arc.loess$reg.type))
arc.type <- arc.type |> 
    arrange(desc(Freq)) |>
    mutate(perc = round(100 * (Freq/sum(Freq)), 2)) 

arc.peaks <- ggplot(arc.type, aes(x="", y=perc, fill=Var1)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0) +
    theme_void() +
    geom_text(aes(label = paste(perc, "%", sep = "")),
              position = position_stack(vjust = 0.5),
              color = "white", size=2) +
    scale_fill_manual(values=cols) +
    labs(fill = "genomic region")

## NT peaks
nt.loess <- peaks.loess[peaks.loess$nt.FDR <= 0.05 & peaks.loess$nt.logFC > 4]
## Format dataframe for plotting
nt.type <- data.frame(table(nt.loess$reg.type))
nt.type <- nt.type |>
    arrange(desc(Freq)) |>
    mutate(perc = round(100 * (Freq/sum(Freq)), 2)) 

nt.peaks <- ggplot(nt.type, aes(x="", y=perc, fill=Var1)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0) +
    theme_void() +
    geom_text(aes(label = paste(perc, "%", sep = "")),
              position = position_stack(vjust = 0.5),
              color = "white", size=2) +
    scale_fill_manual(values=cols) +
    labs(fill = "genomic region")

pdf("./figures/peak_region_types_logFC4.pdf", width = 5, height = 3)
all.peaks + labs(title = paste("all", sum(all.type$Freq), "peaks", sep = " "))
arc.peaks + labs(title = paste(sum(arc.type$Freq), "Arc peaks", sep = " "))
nt.peaks + labs(title = paste(sum(nt.type$Freq), "NT peaks", sep = " "))
dev.off()

