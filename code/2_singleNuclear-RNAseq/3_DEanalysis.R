## #############################################################################
## Date:        February 2025
## Author:      Allison M. Burns
## Filename:    3_DEanalysis.R
## Project:     Cas9 Parse Sequencing
## Description: Run DE analysis on DG excitatory cells using Seurat's
##              FindMakers(). Then analyze gene pathways that are enriched in
##              up-regulated genes
## #############################################################################

################################################################################
## Setup
################################################################################
library(Seurat)
library(tidyverse)
library(gprofiler2)

## Load seurat object
seu <- readRDS("./data/seurat_analysis/SeuratObject_celltype_cas9.rds")
seu <- subset(seu, subclass_name == "037 DG Glut") ## Only get DG glut cells
seu <- subset(seu, sample != "test_1k_") ## remove test samples

## Set parameters
pval = 0.01
fc = 0
org = "mmusculus"

################################################################################
## Run Seurat DE analysis
################################################################################
## Get Cas9+ cells for each sample
seu.cas <- subset(seu, cas9_counts > 0)
Idents(seu.cas) <- seu.cas$orig.ident

## Run DE
de <- FindMarkers(seu.cas,
                  group.by = "orig.ident",
                  slot = "data",
                  ident.1 = "arc",
                  ident.2 = "nt",
                  min.pct = 0.01,
                  logfc.threshold = 0,
                  test.use = "LR",
                  latent.var = "percent.mt")

## Add gene names
de <- rownames_to_column(de, var = "gene")

## Define colors for plotting
de$is.de <- "not de"
de[de$p_val_adj <= pval & de$avg_log2FC < fc,]$is.de <- "down-regulated"
de[de$p_val_adj <= pval & de$avg_log2FC > fc,]$is.de <- "up-regulated"

## Create Volcano Plot
ggplot(de, aes(x = avg_log2FC, y = -log10(p_val), color = is.de)) +
    geom_point() +
    theme_bw()

## Format table
gn2desc <- gconvert(de$gene, organism = org)
de$gn_description <- NA
de$gn_description[gn2desc$input_number] <- gn2desc$description

de <- mutate(de,
             p_val = round(p_val, 2),
             avg_log2FC = round(avg_log2FC,2)) 

################################################################################
## Ontology Analysis
################################################################################
## Run Gene Ontology term enrichment analysis
getGost <- function(x, direction) {
    ## Define list of genes for over-representation analysis
    res <- na.omit(x)
    if(direction == "up") {
        query <- res[res$avg_log2FC > fc & res$p_val_adj <= pval,]
    } else if (direction == "dw") {
        query <- res[res$avg_log2FC < -fc & res$p_val_adj <= pval,]
    }
    query <- query[order(query$p_val_adj, decreasing = FALSE),] ## order genes by p-value
    query <- rownames(query)
    query <- na.omit(query)

    if(length(query) == 0){ query <- "none"}
    
    ## Define list of genes for comparison
    universe <- na.omit(rownames(de))
    
    ## Run Gene Enrichment Analysis
    myGost <- gost(query = query,             
                   organism = org,
                   ordered_query = TRUE,      
                   significant = TRUE,        
                   user_threshold = 0.05,     
                   correction_method = "fdr", 
                   sources = c("GO","KEGG"),  
                   evcodes = TRUE,
                   domain_scope = "custom",   
                   custom_bg = universe)      
    myGost
}

getGostTable <- function(x) {
    ## reformat tables for excel visualization
    if(is.null(x$result)) {
        data.frame(source = as.character(),
                   p_value = as.numeric(),
                   term_name = as.character(),
                   term_size = as.numeric(),
                   query_size = as.numeric(),
                   intersection_size = as.numeric(),
                   intersection = as.character())
    }else{
        x$result |>
            select(source, p_value, term_name, term_size,
                   query_size, intersection_size, intersection) |>
            filter(p_value <= 0.05) |>
            mutate(p_value = round(p_value, 3)) 
    }
}

## Run ontology analysis for up-regulated genes
up.gost <- getGost(de, direction = "up")
up.table <- getGostTable(up.gost)

#################################################################################
## Save tables
#################################################################################
## Save DE gene lists
writexl::write_xlsx(de, path = "./data/DEgenes_seurat.xlsx") 
## Save pathways of up-regulated genes
writexl::write_xlsx(up.table, path = "./data/DEgenes_upreg_KEGG_GO.xlsx") 
