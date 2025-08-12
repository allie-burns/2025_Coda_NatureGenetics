## #############################################################################
## Date:        February 2025
## Author:      Allison M. Burns
## Filename:    2_seurat_processing.R
## Project:     Cas9 Parse Sequencing
## Description: Create and filter seurat object, define cell types for each cell
##              and characterize Cas9 and EGFP positive cells
## #############################################################################

################################################################################
## Setup
################################################################################
library(tidyverse) ## version 2.0.0
library(Seurat) ## version 5.2.1
library(patchwork) ## version 1.3.0

## Set parameters
set.seed(775)
options(future.globals.maxSize= 1891289600)

## Read in the count matrix
mat_path <- "./data/alignment/sublib_comb/all-sample/DGE_filtered/"
mat <- ReadParseBio(mat_path)

## Read in meta data
cell_meta <- read.csv("./data/alignment/sublib_comb/all-sample/DGE_filtered/cell_metadata.csv",
                      row.names = 1)

################################################################################
## Create Seurat Object
################################################################################
## Create seurat object
seu <- CreateSeuratObject(mat,
                          min.features = 100,
                          min.cells = 20,
                          names.field = 0,
                          meta.data = cell_meta)

## Add some information about samples to meta table
seu[["orig.ident"]] <- gsub("[0-9]+_", "", seu$sample)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")

## Visualize raw cell stats
VlnPlot(object   = seu,
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
        group.by = "sample",
        pt.size  = 0,
        ncol     = 3)

wrap_plots(
    FeatureScatter(seu, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5) +
    NoLegend(),
    FeatureScatter(seu, "nCount_RNA","percent.mt", group.by = "orig.ident", pt.size = 0.5) +
    NoLegend(),
    FeatureScatter(seu, "nFeature_RNA","percent.mt", group.by = "orig.ident", pt.size = 0.5),
    ncol = 3) +
    plot_annotation(
        title = "Correlations of cell characteristics - before filtering",
        )

################################################################################
##  Filter Seurat Object
###############################################################################
## Set filter thresholds
min.UMI  <-  500
max.UMI <- 100000
min.genes <- 300
max.genes <- 10000 
max.mito <- 1

## Filter cells by numbers of genes and % mitochondria
seu <- subset(x = seu,
              subset = nFeature_RNA > min.genes & nFeature_RNA < max.genes &
                  nCount_RNA < max.UMI & nCount_RNA > min.UMI &
                  percent.mt < max.mito)

## Visualize filtered cell stats
VlnPlot(object   = seu,
        features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
        group.by = "sample", 
        pt.size  = 0,
        ncol     = 3)

wrap_plots(
    FeatureScatter(seu, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5) +
    NoLegend(),
    FeatureScatter(seu, "nCount_RNA","percent.mt", group.by = "orig.ident", pt.size = 0.5) +
    NoLegend(),
    FeatureScatter(seu, "nFeature_RNA","percent.mt", group.by = "orig.ident", pt.size = 0.5),
    ncol = 3) +
    plot_annotation(
        title = "Correlations of cell characteristics - after filtering",
        )

################################################################################
##  Normalize and integrate Seurat Object
###############################################################################
## Split Seu object into H2 and H28
seu.list <- SplitObject(seu, split.by = "sample")

## normalize and identify variable features for each dataset independently
seu.list <- lapply(X = seu.list, FUN = function(x) {
    ## Normalize data with sctransform
    x <- SCTransform(object = x,
                     method = "glmGamPoi",
                     vars.to.regress = "percent.mt",
                     variable.features.n = 3000)
    ## Normalize data with log normalization
    DefaultAssay(x) <- "RNA"
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
    ## Save object for integration
    DefaultAssay(x) <- "SCT"
    x
})

## perform integration
features <- SelectIntegrationFeatures(object.list = seu.list)
seu.anchors <- FindIntegrationAnchors(object.list = seu.list, anchor.features = features)
seu.combined <- IntegrateData(anchorset = seu.anchors)

################################################################################
## Dimensionality reduction and Clustering
################################################################################
## Cluster based on integrated data
DefaultAssay(seu.combined) <- "integrated"

## Run the standard workflow for visualization and clustering
seu.combined <- ScaleData(seu.combined, verbose = FALSE)
seu.combined <- RunPCA(seu.combined, npcs = 50, verbose = FALSE)

## Decide on PCS
ElbowPlot(seu.combined, ndims = 50)

## Run TSNE and UMAP and find clusters
seu.combined <- RunTSNE(seu.combined, dims = 1:20, verbose = FALSE)
seu.combined <- RunUMAP(seu.combined, reduction = "pca", dims = 1:20)
seu.combined <- FindNeighbors(seu.combined, reduction = "pca", dims = 1:20)
seu.combined <- FindClusters(seu.combined, resolution = 0.8)

## Visualization
DimPlot(seu.combined, group.by = "orig.ident", shuffle = TRUE)
DimPlot(seu.combined, group.by = "sample", shuffle = TRUE)
DimPlot(seu.combined, label = TRUE)

################################################################################
## Check expression on UMAP
################################################################################
FeaturePlot(seu.combined,
            features = c("nCount_RNA","nFeature_RNA","percent.mt"),
            order = TRUE,
            ncol = 3)

FeaturePlot(seu.combined,
            features = c("nCount_SCT","nFeature_SCT","percent.mt"),
            order = TRUE,
            ncol = 3)

################################################################################
## Annotate CellTypes with MapMyCells
################################################################################
## Read in Seurat object
seu <- seu.combined

## Add count data to RNA assay
DefaultAssay(seu) <- "RNA"
seu@assays$RNA$counts <- seu@assays$SCT$counts
seu <- NormalizeData(seu)

## Count cells in each seurat cluster
cellCount <- table(seu@meta.data$seurat_clusters, seu@meta.data$orig.ident)
cellCount <- data.frame(cellCount)
cellCount |>
    spread(key = Var2, value = Freq) |>
    mutate(total = arc + nt) |>
    mutate(perc = paste(round(100*(total/sum(total)),1),"%",sep = ""))

## prepare data for MapMyCells
rnaCounts <- as.matrix(seu@assays$RNA$counts) ## Read count matrix (raw gene activities)
obs = colnames(rnaCounts) # Store sample-based metadata in a dataframe called obs
var2 = rownames(rnaCounts) # gene data 
obs <- data.frame(obs)
rownames(obs) <- obs$obs
var2 <- data.frame(var2)
rownames(var2) <- var2$var2
count_matrix <- t(rnaCounts) ## Transpose counts so samples are rows and genes are columns
count_matrix <- as.sparse(count_matrix)
ad <- anndata::AnnData(X = count_matrix, obs = obs, var = var2)
output_path = "./data/adata_for_MapMyCells.h5ad" # Write with compression
anndata::write_h5ad(ad, output_path, compression='gzip')
file_size_bytes = file.size(output_path) # Determine and print the file size
print(paste("File size in bytes:", file_size_bytes))

## From here, load file into MapMyCells (https://portal.brain-map.org/atlases-and-data/bkp/mapmycells) and then download the output. This file can then be loaded in and added to the Seurat object

## merge mapmycells output with seurat object
mapmycells <- read.csv("./data/seurat_analysis/adata_for_MapMyCells_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1739960896951/adata_for_MapMyCells_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1739960896951.csv", skip = 4)
## Merge data.frames
mmc <- select(mapmycells,
              c("cell_id",
                "class_name","class_bootstrapping_probability",
                "subclass_name", "subclass_bootstrapping_probability",
                "supertype_name","supertype_bootstrapping_probability",
                "cluster_name","cluster_bootstrapping_probability"))
seu@meta.data <- cbind(seu@meta.data,mmc[match(mmc$cell_id, Cells(seu)),])

## Visualize new cell type assignments
FeaturePlot(seu,
            features = c("class_bootstrapping_probability",
                         "subclass_bootstrapping_probability",
                         "supertype_bootstrapping_probability",
                         "cluster_bootstrapping_probability"),
            order = FALSE)

## Visualize Cell types
DimPlot(seu, group.by = "class_name")
DimPlot(seu, group.by = "class_name",label = TRUE)
DimPlot(seu, group.by = "subclass_name")
DimPlot(seu, group.by = "subclass_name",label = TRUE)
DimPlot(seu, group.by = "supertype_name")
DimPlot(seu, group.by = "supertype_name",label = TRUE)
DimPlot(seu, group.by = "cluster_name")
DimPlot(seu, group.by = "cluster_name",label = TRUE)

## What happens if I remove low probability and cell types with few counts for subclass?
subclass.count <- data.frame(table(seu$subclass_name))
seu.subclass <- subset(seu,
                       subclass_bootstrapping_probability == 1 &
                       subclass_name %in% subclass.count$Var1[subclass.count$Freq > 10])
DimPlot(seu.subclass, group.by = "subclass_name")
DimPlot(seu.subclass, group.by = "subclass_name",label = TRUE)

## Subset meta data object to only include relevant columns
seu.subclass@meta.data <- select(
    seu.subclass@meta.data,
    c("orig.ident", "sample", "nCount_RNA", "nFeature_RNA", "percent.mt",
      "nCount_SCT", "nFeature_SCT", "seurat_clusters", "cell_id",
      "subclass_name","subclass_bootstrapping_probability"))

################################################################################
## Define Cas9 and EGFP positive cells
################################################################################
## Annotated cells as seurat object
seu <- seu.subclass

## Load counts table
mat_path <- "./data/alignment/sublib_comb/all-sample/DGE_filtered/"
mat <- ReadParseBio(mat_path)

## Add Cas9 and EGFP raw counts to meta.data
seu.counts <- mat[,match(seu$cell_id,colnames(mat))]
seu$cas9_counts <- seu.counts[rownames(seu.counts) == "Cas9m4",]
seu$egfp_counts <- seu.counts[rownames(seu.counts) == "EGFP",]

## Add Cas9 raw counts to counts table
tab <- rbind(seu.counts[rownames(seu.counts) == "Cas9m4",],seu@assays$RNA$counts)
rownames(tab)[rownames(tab) == ""] <-  "Cas9m4"
seu@assays$RNA$counts <- tab
seu <- NormalizeData(seu) ## Add new normalized value to data

## Plot Cas9 and EGFP + cells
FeaturePlot(seu,
            features = c("cas9_counts","egfp_counts"),
            order = TRUE)

## Count cas9 and EGFP+cells
pos_counts <- table(seu$subclass_name, seu$cas9_counts > 0, seu$orig.ident) |>
    data.frame() |>
    filter(Var2 == "TRUE") |>
    group_by(Var1) |>
    filter(sum(Freq) != 0)

## Number of Cas9+ cells in each cell type
ggplot(pos_counts, aes(x = Var1, y = Freq, fill = Var3)) +
    geom_bar(stat = "identity", position= position_dodge()) +
    geom_text(aes(label = Freq), position = position_dodge(width = 0.9), vjust = -0.2) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Compare raw and normalized counts
VlnPlot(seu, features = c("Cas9m4","EGFP"), group.by = "sample", assay = "RNA",slot = "counts")
VlnPlot(seu, features = c("Cas9m4","EGFP"), group.by = "sample", assay = "RNA",slot = "data")

## Save seurat object with cell type annotations and Cas9 and EGFP counts defined
saveRDS(seu,"./data/seurat_analysis/SeuratObject_celltype_cas9.rds")
