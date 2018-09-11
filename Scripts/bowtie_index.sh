#!/bin/bash

#***************************************************************#
#         Generate the Salmon Genome Index using Bowtie         #
#                 STEP 1 IN BOWTIE-TOPHAT PIPELINE              #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=BOWTIE_index
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Logs/bowtie_%A.out

module load bowtie/2.2.9

#Working Directories
INPDIR=/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Genomes/Salmo_salar
OUTDIR=/flush3/sam079/RNAseq-POMV/GenomeIndex/Bowtie

bowtie2-build --threads 8 ${INPDIR}/GCF_000233375.1_ICSASG_v2_genomic.fna ${OUTDIR}/Salmon
