#!/bin/bash

#***************************************************************#
#              Assemble Transcripts with Cufflinks              #
#           STEP 3 IN BOWTIE-TOPHAT-CUFFLINKS PIPELINE          #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=CUFFLINKS
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Logs/cufflinks_%A.out


module load cufflinks/2.2.2

# Working Directories

INPDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment/AlignSalmonGenomeBowtie
OUTDIR=/flush3/sam079/RNAseq-POMV/Processed

cufflinks -p 8 -o ${OUTDIR}/TEST ${INPDIR}/POMV6HPIR1/accepted_hits.bam


