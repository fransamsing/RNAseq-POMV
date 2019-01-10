#!/bin/bash

#***************************************************************#
#            Generate Viral Genome Index POMV and ISA           #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=bowtie_index
#SBATCH --time=00:00:10
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Logs/bowtie_index_%A.out


module load bowtie/2.2.9

# Working Dictories 
INPDIR=/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Genomes/POMV_genome_AAHL
#INPDIR=/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Genomes/ISA_genome

#OUTDIR=/flush3/sam079/RNAseq-POMV/ViralGenomeIndex/POMVGenomeIndex
#OUTDIR=/flush3/sam079/RNAseq-POMV/ViralGenomeIndex/ISAVGenomeIndex
OUTDIR=/flush3/sam079/RNAseq-POMV/ViralGenomeIndex/POMVNewGenomeIndex


#bowtie2-build ${INPDIR}/POMV_14_01514_ORF.fa ${OUTDIR}/POMV
#bowtie2-build ${INPDIR}/ISA_Glesvaer_2_90_ORF.fa ${OUTDIR}/ISAV
bowtie2-build ${INPDIR}/POMV_14_01514_GENOME_format.fa ${OUTDIR}/POMV