#!/bin/bash

#***************************************************************#
#              Assemble Transcripts with Cufflinks              #
#           STEP 3 IN BOWTIE-TOPHAT-CUFFLINKS PIPELINE          #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=CUFFLINKS
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=15GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Logs/cufflinks_%A.out
#SBATCH --array=0-14

module load jemalloc
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE
module load cufflinks

# Working Directories

INPDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment/AlignSalmonPOMVCombined
OUTDIR=/flush3/sam079/RNAseq-POMV/Processed/Assembly/SalmonPOMV

SAMPLES=( $(cut -d , -f 1 ../STARInputList.csv) );


if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
then
    i=$SLURM_ARRAY_TASK_ID
    INFILES=${INPDIR}/${SAMPLES[$i]}
    cufflinks -p 8 -o ${SAMPLES[$i]} ${INFILES}/accepted_hits.bam
else
    echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi

# I could not create a directory for the output files so I just run 
# this script and then manually moved the files 
# to the flush directory in OUTDIR 


