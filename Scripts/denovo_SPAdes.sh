#!/bin/bash

#***************************************************************#
#                De novo assembly of POMV using spades   	    #                 
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=spades_denovo_assembly
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Logs/spades_%A.out

module load spades/3.12.0

INPDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment/AlignSalmonGenome/Star
OUTDIR=/flush3/sam079/RNAseq-POMV/Processed/DenovoAssemblyPOMV

#cp ${INPDIR}/POMV24HPIR1Unmapped.out.mate1 ${INPDIR}/POMV24HPIR1Unmapped.out.mate1.fastq
#cp ${INPDIR}/POMV24HPIR1Unmapped.out.mate2 ${INPDIR}/POMV24HPIR1Unmapped.out.mate2.fastq

spades.py --rna -1 ${INPDIR}/POMV24HPIR1Unmapped.out.mate1.fastq -2 ${INPDIR}/POMV24HPIR1Unmapped.out.mate2.fastq -o ${OUTDIR}/SPAdes

