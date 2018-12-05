#!/bin/bash

#***************************************************************#
#              Using samtools to get alignment stats            #
#               		for viral genomes  			            #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=samtools_sort
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Logs/samtools_sort_%A.out
#SBATCH --array=0-14

#-----------------------------------------------------------------#

module load jemalloc
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE
module load samtools/0.1.19
module load bedtools/2.26.0

# Working Directories
INPDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment/AlignPOMVGenome
SAMPLES=( $(cut -d , -f 1 ../STARInputList.csv) );

if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
then
    i=$SLURM_ARRAY_TASK_ID
	samtools view -bS ${INPDIR}/${SAMPLES[$i]}.sam > ${INPDIR}/${SAMPLES[$i]}.bam
	samtools sort ${INPDIR}/${SAMPLES[$i]}.bam ${INPDIR}/${SAMPLES[$i]}.sorted
	bamToBed -i ${INPDIR}/${SAMPLES[$i]}.sorted.bam > ${INPDIR}/${SAMPLES[$i]}.sorted.bed
	genomeCoverageBed -i ${INPDIR}/${SAMPLES[$i]}.sorted.bed -g ../POMV_sequences.txt -dz > ${INPDIR}/${SAMPLES[$i]}_coverage.txt
else
    echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi