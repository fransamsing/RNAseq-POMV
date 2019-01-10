#!/bin/bash

#***************************************************************#
#                De novo assembly of POMV using TRINITY	 	    #                 
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=Trinity_denovo_assembly
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Logs/Trinity_%A.out

module load trinity/2.8.4

INPDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment/AlignSalmonGenome/Star
OUTDIR=/flush3/sam079/RNAseq-POMV/Processed/DenovoAssemblyPOMV

## To run Trinity, reads have to be formated in the old Ilumina format. 
## STAR does not have an option to add these marks, however, you can easily add them - separately for each of the 
## Unmapped mate file using:
# awk '{if (NR%4==1) $1=$1 "/1"; print}' ${INPDIR}/POMV24HPIR1Unmapped.out.mate1.fastq > ${INPDIR}/Unmapped.out.mate1.mark12
# awk '{if (NR%4==1) $1=$1 "/2"; print}' ${INPDIR}/POMV24HPIR1Unmapped.out.mate2.fastq > ${INPDIR}/Unmapped.out.mate2.mark12

Trinity --seqType fq  --SS_lib_type RF --max_memory 40G --min_kmer_cov 1 --CPU 8 \
--left ${INPDIR}/Unmapped.out.mate1.mark12 \
--right ${INPDIR}/Unmapped.out.mate2.mark12 \
--output ${OUTDIR}/Trinity

## .fastq files were created using the following. This is done to run SPADES.
## cp ${INPDIR}/POMV24HPIR1Unmapped.out.mate1 ${INPDIR}/POMV24HPIR1Unmapped.out.mate1.fastq
## cp ${INPDIR}/POMV24HPIR1Unmapped.out.mate2 ${INPDIR}/POMV24HPIR1Unmapped.out.mate2.fastq







