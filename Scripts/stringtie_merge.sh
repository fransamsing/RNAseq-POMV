#!/bin/bash

#***************************************************************#
#              Assemble Transcripts with StringTie              #
#           STEP 3 IN BOWTIE-TOPHAT-STRINGTIE PIPELINE          #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=STRINGTIE
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=15GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Logs/stringtie_%A.out

module load stringtie

# Working Directories
ANODIR=/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Genomes/Salmo_salar
INPDIR=/flush3/sam079/RNAseq-POMV/Processed/Assembly/SalmonPOMV/StringTie
OUTDIR=/flush3/sam079/RNAseq-POMV/Processed/Assembly/SalmonPOMV


stringtie --merge -p 8 -G ${ANODIR}/GCF_000233375.1_ICSASG_v2_genomic.gff -o ${OUTDIR}/SalmonPOMV_stringtie_merged.gtf ../all_SalmonPOMV_transcripts.txt


