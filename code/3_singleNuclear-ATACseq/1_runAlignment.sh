#!/bin/bash

## #############################################################################
## Date:        October 2024
## Author:      Allison M. Burns
## Filename:    1_runAlignment.sh
## Project:     single cell ATAC seq
## Description: Run alignments using cellranger
## #############################################################################

################################################################################
## Align reads
################################################################################
cellranger-atac count \
		--id=dCas9_plus_Arc_sgRNA \
		--reference=./CRREF025_ARC-mm10_with_Cas9m4_EGFP/mm10 \
		--fastqs=./fastq_files/arc_ \
		--sample=dCas9_plus_Arc_sgRNA

cellranger-atac count \
		--id=dCas9_plus_NT_sgRNA \
		--reference=./CRREF025_ARC-mm10_with_Cas9m4_EGFP/mm10 \
		--fastqs=./fastq_files/nt_ \
		--sample=dCas9_plus_NT_sgRNA

## Subset reads aligning to EGFP or Cas9 from bam files
samtools view -b ./cellranger/dCas9_plus_Arc_sgRNA/outs/possorted_bam.bam EGFP > ./subset_bam/dCas9_Arc_egfp.bam
samtools view -b ./cellranger/dCas9_plus_Arc_sgRNA/outs/possorted_bam.bam Cas9m4 > ./subset_bam/dCas9_Arc_cas9m4.bam

samtools view -b ./cellranger/dCas9_plus_NT_sgRNA/outs/possorted_bam.bam EGFP > ./subset_bam/dCas9_NT_egfp.bam
samtools view -b ./cellranger/dCas9_plus_NT_sgRNA/outs/possorted_bam.bam Cas9m4 > ./subset_bam/dCas9_NT_cas9m4.bam

