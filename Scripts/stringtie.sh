#!/bin/bash

#***************************************************************#
#              Assemble Transcripts with StringTie              #
#           STEP 3 IN BOWTIE-TOPHAT-STRINGTIE PIPELINE          #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=STRINGTIE
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=15GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Logs/stringtie_%A.out
#SBATCH --array=0-14

module load jemalloc
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE
module load stringtie

# Working Directories

INPDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment/AlignSalmonPOMVCombined
OUTDIR=/flush3/sam079/RNAseq-POMV/Processed/Assembly/SalmonPOMV/StringTie

SAMPLES=( $(cut -d , -f 1 ../STARInputList.csv) );


if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
then
    i=$SLURM_ARRAY_TASK_ID
    INFILES=${INPDIR}/${SAMPLES[$i]}
    OUTFILES=${OUTDIR}/${SAMPLES[$i]}
    stringtie ${INFILES}/accepted_hits.bam -p 8 -o ${OUTFILES} 
else
    echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi
