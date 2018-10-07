#!/bin/bash

#***************************************************************#
#                  Building read counts matrix                  #
#             STEP XX in STAR-FeatureCounts pipeline            #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=FeactureCounts
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Logs/feature_counts_%A.out

module load subread

# Working Directories
REFDIR=/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Genomes/Salmo_salar 
INPDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment/AlignSalmonGenomeStar/STAR
OUTDIR=/flush3/sam079/RNAseq-POMV/Processed/Expression/Read_counts

featureCounts -T 8 -p -t exon \
-g gene \
-a ${REFDIR}/GCF_000233375.1_ICSASG_v2_genomic.gff \
-o ${OUTDIR}/read_counts.txt \
${INPDIR}/*Aligned.sortedByCoord.out.bam \ 






