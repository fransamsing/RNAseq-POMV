#!/bin/bash

#***************************************************************#
#       Align denovo assembled viruses to available genome      #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=bowtie_align_denovo
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=1GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Logs/bowtie_align_%A.out

#---------------------------------------------------------------#

module load gmap/2017-11-15 
module load samtools/1.3.1
module load bowtie

REFDIR=/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Genomes/POMV_genome_AAHL 
READSDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment/AlignSalmonGenome/Star
OUTDIR=/flush3/sam079/RNAseq-POMV/Processed/DenovoAssemblyPOMV

## POMV ORFs
IDXDIR=/flush3/sam079/RNAseq-POMV/ViralGenomeIndex/POMVGenomeIndex

bowtie2 -p 8 -x ${IDXDIR}/POMV \
-1 ${READSDIR}/POMV24HPIR1Unmapped.out.mate1.fastq \
-2 ${READSDIR}/POMV24HPIR1Unmapped.out.mate2.fastq  \
-S ${OUTDIR}/BOWTIE_unmapped.sam 2> ${OUTDIR}/BOWTIE_unmapped.log 

#REFDIR=/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Genomes/ISA_genome
#READSDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment/AlignSalmonGenome/Star
#OUTDIR=/flush3/sam079/RNAseq-POMV/Processed/DenovoAssemblyISAV

## ISAV ORFs
#IDXDIR=/flush3/sam079/RNAseq-POMV/ViralGenomeIndex/ISAVORFIndex 

#bowtie2 -p 8 -x ${IDXDIR}/ISAV \
#-1 ${READSDIR}/ISAV24HPIR1Unmapped.out.mate1.fastq \
#-2 ${READSDIR}/ISAV24HPIR1Unmapped.out.mate2.fastq  \
#-S ${OUTDIR}/BOWTIE_unmapped.sam 2> ${OUTDIR}/BOWTIE_unmapped.log 

samtools view -Sb ${OUTDIR}/BOWTIE_unmapped.sam > ${OUTDIR}/BOWTIE_unmapped.bam \

samtools sort ${OUTDIR}/BOWTIE_unmapped.bam -o ${OUTDIR}/BOWTIE_unmapped.sorted.bam \

samtools index ${OUTDIR}/BOWTIE_unmapped.sorted.bam

## POMV New Genome
#IDXDIR=/flush3/sam079/RNAseq-POMV/ViralGenomeIndex/POMVNewGenomeIndex

#bowtie2 -p 8 -x ${IDXDIR}/POMV \
#-1 ${READSDIR}/POMV24HPIR1Unmapped.out.mate1.fastq \
#-2 ${READSDIR}/POMV24HPIR1Unmapped.out.mate2.fastq  \
#-S ${OUTDIR}/BOWTIE_new_unmapped.sam 2> ${OUTDIR}/BOWTIE_new_unmapped.log 

#samtools view -Sb ${OUTDIR}/BOWTIE_new_unmapped.sam > ${OUTDIR}/BOWTIE_new_unmapped.bam \

#samtools sort ${OUTDIR}/BOWTIE_new_unmapped.bam -o ${OUTDIR}/BOWTIE_new_unmapped.sorted.bam \

#samtools index ${OUTDIR}/BOWTIE_new_unmapped.sorted.bam 