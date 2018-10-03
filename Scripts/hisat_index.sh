#!/bin/bash

#***************************************************************#
#         Generate the salmon Genome Index using HISAT          #
#                 STEP 1 IN HISAT PIPELINE                      #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=HISAT_index
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Logs/HISAT_index_%A.out

module load hisat/2.0.5

#Working Directories
INPDIR=/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Genomes
OUTDIR=/flush3/sam079/RNAseq-POMV/GenomeIndex/Hisat_SalmonPOMV

#cat $INPDIR/Salmo_salar/GCF_000233375.1_ICSASG_v2_genomic.fna $INPDIR/POMV_genome_AAHL/POMV_14_01514_ORF.fa > $INPDIR/SalmonPOMV_genomic.fa

hisat2-build --threads 8 ${INPDIR}/SalmonPOMV_genomic.fa ${OUTDIR}/SalmonPOMV 



