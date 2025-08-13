# 2025_Coda_NatureGenetics

https://doi.org/10.5281/zenodo.16834060

This is the code used to perform the analyses in the manuscript, "Cell-type and locus-specific epigenetic editing of memory expression" (Coda et al., 2025, Nature Genetics).


## Abstract
Epigenetic mechanisms have for long been proposed to act as molecular mnemonics1–3, but whether the epigenetic makeup of a single genomic site can guide learnt behaviors remains unknown. Here, we combined CRISPR-based epigenetic editing tools4,5 with c-Fos driven engram technologies6,7 to address this question in memory-bearing neuronal ensembles. Focusing on the promoter region of Arc, a master regulator of synaptic plasticity8, we found that its locus-specific and temporally controllable epigenetic editing is necessary and sufficient to regulate memory expression. Such effects occurred irrespective of the memory phase – during the initially labile period after learning and for fully consolidated memories – and were reversible within subject, testifying to their inherent plasticity. These findings provide a proof-of-principle that site-specific epigenetic dynamics are causally implicated in memory expression.


## Data Availability
Raw and processed data files can be downloaded from [GSE299742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299742).


## Analysis Information
### Bulk ChIP-sequencing
This analysis was performed in order to determine off-target effects of Cas9.

Fastq files ([GSE256419](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE256419)) were aligned and processed using Bowtie2, samtools and picard (`/code/1_bulk_ChIPseq/1_runAlignment.sh`). The R tool, csaw, was used to call peaks (`/code/1_bulk_ChIPseq/2_PeakCalling_csaw.R`).  And Cas9 targets were determined by performing a differential enrichment analysis between IP and input samples (`/code/1_bulk_ChIPseq/3_DEanalysis_edgeR.R`)


### Single-cell RNA-sequencing
This analysis was performed to see which genes and pathways are enriched in Cas9+ excitatory neurons of the DG after Arc inactivation.

Fastq files ([GSE299740](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299740)) were aligned and processed using Parse's Trailmaker (`/code/2_singleNuclear-RNAseq/1_runAlignment.sh`). Seurat was used to process the counts data, define cell types and characterize the cells (`/code/2_singleNuclear-RNAseq/2_seurat_processing.R`). Finally, Seurat was used to define genes that were differentially expressed in the Cas9 positive cells of the Dentate Gyrus (`/code/2_singleNuclear-RNAseq/3_DEanalysis.R`). 


### Single-cell ATAC-sequencing
This analysis was performed to see which regions have differential expression after Cas9 targeting of the Arc promoter. 

Fastq files ([GSE299741](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299741)) were aligned and processed using CellRanger (`/code/3_singleNuclear-ATACseq/1_runAlignment.sh`). Seurat and Signac were used to process the counts data, define cell types, characterize the cells and perform differential expression analysis (`/code/3_singleNuclear-ATACseq/2_SignacAnalysis.R`). 
