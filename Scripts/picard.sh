#!/bin/bash

#***************************************************************#
#            Picard to extract Alignment Summary Metrics        # 
#                    Runs Picard on BAM files                   #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=picard
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Logs/picard_%A.out

module load picard/2.18.2 

java -jar -Xmx5g /apps/picard/2.18.2/picard.jar CollectAlignmentSummaryMetrics \
METRIC_ACCUMULATION_LEVEL=ALL_READS \
R=/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Genomes/Salmo_salar/GCF_000233375.1_ICSASG_v2_genomic.fna \
INPUT=/flush3/sam079/RNAseq-POMV/Processed/Alignment/AlignSalmonGenomeBowtie/TophatDefaults/NEGcontrolR1/accepted_hits.bam \
OUTPUT=/flush3/sam079/RNAseq-POMV/Processed/output_picard_test_tophat.txt 
