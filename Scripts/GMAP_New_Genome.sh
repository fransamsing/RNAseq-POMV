#!/bin/bash

#***************************************************************#
#       Align denovo assembled contigs using GMAP               #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=GMAP
#SBATCH --time=00:02:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Logs/GMAP_%A.out

#---------------------------------------------------------------#
module load gmap/2017-11-15 
module load samtools/1.3.1

REFDIR=/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Genomes/POMV_genome_AAHL 
INPDIR=/flush3/sam079/RNAseq-POMV/Processed/DenovoAssemblyPOMV
READSDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment/AlignSalmonGenome/Star
OUTDIR=/flush3/sam079/RNAseq-POMV/Processed/DenovoAssemblyPOMV

## GMAP NEW INDEX
gmap_build -d POMV_new_genome -D ${OUTDIR}/GMAP \
-k 13 ${REFDIR}/POMV_14_01514_New_GENOME_format.fa

## SPADES - NEW POMV GENOME
gmap -n 0 -D ${OUTDIR}/GMAP \
-d POMV_new_genome ${INPDIR}/SPAdes/transcripts.fasta \
-f samse > ${OUTDIR}/GMAP/spades_new_gmap.sam

samtools view -Sb ${OUTDIR}/GMAP/spades_new_gmap.sam > ${OUTDIR}/GMAP/spades_new_gmap.bam \

samtools sort ${OUTDIR}/GMAP/spades_new_gmap.bam -o ${OUTDIR}/GMAP/spades_new_gmap.bam \

samtools index ${OUTDIR}/GMAP/spades_new_gmap.bam

## TRINITY - NEW POMV GENOME
gmap -n 0 -D ${OUTDIR}/GMAP \
-d POMV_new_genome ${INPDIR}/Trinity/Trinity.fasta \
-f samse > ${OUTDIR}/GMAP/trinity_new_gmap.sam

samtools view -Sb ${OUTDIR}/GMAP/trinity_new_gmap.sam > ${OUTDIR}/GMAP/trinity_new_gmap.bam \

samtools sort ${OUTDIR}/GMAP/trinity_new_gmap.bam -o ${OUTDIR}/GMAP/trinity_new_gmap.bam \

samtools index ${OUTDIR}/GMAP/trinity_new_gmap.bam

