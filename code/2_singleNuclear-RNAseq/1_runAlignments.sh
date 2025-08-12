#!/usr/bin/env bash

## #############################################################################
## Date:        February 2025
## Author:      Allison M. Burns
## Filename:    1_runAlignments.sh
## Project:     Cas9 Parse Sequencing
## Description: Run single nuclear alignment using TrailMaker from Parse Bio-
##              Sciences and merge sublibraries. 
## #############################################################################

################################################################################
## Setup
################################################################################
## activate conda environment to run Parse seq
module load gcc python 
source activate /home/burns/miniconda3/envs/spipe
conda activate spipe

## Point to folder
myDir="./Parse_snSeq"

################################################################################
## TrailMaker Alignment for each sublibrary
################################################################################
## List sublibraries
sublib=$(ls -R $myDir/expdata/fastq_files | grep "fastq.gz" | cut -c 1-11 | uniq)

## Run Alignment
for s in $sublib
do
    ## Merge siblibrary files across runs
    r1=$(find $myDir/expdata/fastq_files | grep $s | grep "fastq.gz"| grep "R1")
    r2=$(find $myDir/expdata/fastq_files | grep $s | grep "fastq.gz"| grep "R2")
    ## Cat fastq files from multiple lanes
    if [[ ! -e $myDir/expdata/fastq_files/${s}_S1_R1_001.fastq.gz ]]; then
	echo Merging: $s
	cat $r1 > $myDir/expdata/fastq_files/${s}_S1_R1_001.fastq.gz
	cat $r2 > $myDir/expdata/fastq_files/${s}_S1_R2_001.fastq.gz
    fi

    ## Run Alignment 
    echo Aligning: $s
    split-pipe \
	--mode all \
	--chemistry v3 \
	--genome_dir ./genomes/mm10_cas9 \
	--fq1 $myDir/expdata/fastq_files/${s}_S1_R1_001.fastq.gz \
	--fq2 $myDir/expdata/fastq_files/${s}_S1_R2_001.fastq.gz \
	--output_dir $myDir/analysis/$s \
	--samp_sltab $myDir/analysis/Parse_Biosciences_Evercode_WT_Sample_Loading_Table_v2.xlsm
done

## Merge sublibraries
split-pipe \
    --mode comb \
    --sublibraries \
    $myDir/analysis/sublibrary1 \
    $myDir/analysis/sublibrary2 \
    $myDir/analysis/sublibrary3 \
    $myDir/analysis/sublibrary4 \
    $myDir/analysis/sublibrary5 \
    $myDir/analysis/sublibrary6 \
    $myDir/analysis/sublibrary7 \
    $myDir/analysis/sublibrary8 \
    --output_dir $myDir/analysis/sublib_comb

