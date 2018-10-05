#!/bin/bash

#***************************************************************#
#              HISAT Alignment with Salmon Genome               #
#               STEP 3 IN HISAT-Stringtie pipeline              #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=HISAT_samtobam
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Logs/HISAT_samtobam_%A.out
#SBATCH --array=0-14

module load jemalloc
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE
module load samtools/1.7.0 

# Working Directories
INPDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment/AlignSalmonPOMVHisat
OUTDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment/AlignSalmonPOMVHisat

SAMPLES=( $(cut -d , -f 1 ../STARInputList.csv) );

if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
then
    i=$SLURM_ARRAY_TASK_ID
    INFILES=${INPDIR}/${SAMPLES[$i]}
    samtools sort -@ 8 -o ${OUTDIR}/${SAMPLES[$i]}.bam ${INFILES}.sam
else
    echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi
