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

#REFDIR=/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Genomes/POMV_genome_AAHL 
#INPDIR=/flush3/sam079/RNAseq-POMV/Processed/DenovoAssemblyPOMV
#OUTDIR=/flush3/sam079/RNAseq-POMV/Processed/DenovoAssemblyPOMV

REFDIR=/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Genomes/ISA_genome 
INPDIR=/flush3/sam079/RNAseq-POMV/Processed/DenovoAssemblyISAV
OUTDIR=/flush3/sam079/RNAseq-POMV/Processed/DenovoAssemblyISAV

## GMAP INDEX
#gmap_build -d POMV_genome -D ${OUTDIR}/GMAP \
gmap_build -d ISAV_orfs -D ${OUTDIR}/GMAP \
-k 13 ${REFDIR}/ISA_Glesvaer_2_90_ORF.fa 

## SPADES ## CHANGE -d BETWEEN DIFFERENT VIRAL GENOMES 
gmap -n 0 -D ${OUTDIR}/GMAP \
-d ISAV_orfs ${INPDIR}/SPAdes/transcripts.fasta \
-f samse > ${OUTDIR}/GMAP/spades_gmap.sam

samtools view -Sb ${OUTDIR}/GMAP/spades_gmap.sam > ${OUTDIR}/GMAP/spades_gmap.bam \

samtools sort ${OUTDIR}/GMAP/spades_gmap.bam -o ${OUTDIR}/GMAP/spades_gmap.bam \

samtools index ${OUTDIR}/GMAP/spades_gmap.bam

## TRINITY ## CHANGE -d BETWEEN DIFFERENT VIRAL GENOMES 
gmap -n 0 -D ${OUTDIR}/GMAP \
-d ISAV_orfs ${INPDIR}/Trinity/Trinity.fasta \
-f samse > ${OUTDIR}/GMAP/trinity_gmap.sam

samtools view -Sb ${OUTDIR}/GMAP/trinity_gmap.sam > ${OUTDIR}/GMAP/trinity_gmap.bam \

samtools sort ${OUTDIR}/GMAP/trinity_gmap.bam -o ${OUTDIR}/GMAP/trinity_gmap.bam \

samtools index ${OUTDIR}/GMAP/trinity_gmap.bam