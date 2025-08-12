#!/bin/bash

## #############################################################################
## Date:        March 2023
## Author:      Allison M. Burns
## Filename:    1_runAlignment.sh
## Project:     Cas9 ChIP
## Description: Run alignments, sort and index bam files, mark duplicates, and 
##              get mapping statistics for downstream visualization.
## #############################################################################

## Get list of files
path="/data_directory/"
## Get sample names (for pulling out both paired files)
samples=$(ls $path/fastq_files/*.fastq.gz | sed 's#.*/##' | rev | cut -c20- | rev | uniq)

for samp in $samples
do
    ## Set up file system for STAR runs
    echo Aligning: $samp 
    r1=$(ls $path/fastq_files/*.fastq.gz | grep $samp | grep "R1")
    r2=$(ls $path/fastq_files/*.fastq.gz | grep $samp | grep "R2")
        
    ## Run Bowtie2
    bowtie2 -p 24 -q \
     	    -x /data_directory/genomes/m_musculus/mm10/mm10_bowtie2_build/GrCm38 \
	    -1 $r1\
	    -2 $r2 \
	    -S $path/1_bowtie_align/${samp}_unsorted.sam
    
    ## SAM2BAM
    samtools view -h -S -b \
    	     -o $path/1_bowtie_align/${samp}_unsorted.bam \
    	     $path/1_bowtie_align/${samp}_unsorted.sam
    
    ## sort bam
    samtools sort $path/1_bowtie_align/${samp}_unsorted.bam \
	     -o $path/1_bowtie_align/${samp}_sorted.bam 
    
    ## mark duplicates
    picard MarkDuplicates \
	   I=$path/1_bowtie_align/${samp}_sorted.bam \
	   O=$path/2_mark_dups/${samp}_sorted_mdups.bam \
	   M=$path/2_mark_dups/${samp}_sorted_mdups_metrics.txt

    ## samtools flagstat
    samtools flagstat $path/2_mark_dups/${samp}_sorted_mdups.bam > $path/2_mark_dups/${samp}_flagstat.txt 
done
