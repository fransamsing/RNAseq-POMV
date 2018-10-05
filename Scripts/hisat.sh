#!/bin/bash

#***************************************************************#
#              HISAT Alignment with Salmon Genome               #
#               STEP 2 IN HISAT-Stringtie pipeline              #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=HISAT
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Logs/HISAT_%A.out
#SBATCH --array=0-14

module load jemalloc
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE
module load hisat/2.0.5

# Working Directories
INPDIR=/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Data
REFDIR=/flush3/sam079/RNAseq-POMV/GenomeIndex/Hisat_SalmonPOMV
OUTDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment/AlignSalmonPOMVHisat

SAMPLES=( $(cut -d , -f 1 ../STARInputList.csv) );
INFILES_R1_LIST=( $(cut -d , -f 2 ../STARInputList.csv) );
INFILES_R2_LIST=( $(cut -d , -f 3 ../STARInputList.csv) );

if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
then
    i=$SLURM_ARRAY_TASK_ID
    INFILES_R1=${INPDIR}/${INFILES_R1_LIST[$i]}
    INFILES_R2=${INPDIR}/${INFILES_R2_LIST[$i]}
    hisat2 -p 8 --dta -x ${REFDIR}/SalmonPOMV -S ${OUTDIR}/${SAMPLES[$i]}.sam -1 ${INFILES_R1} -2 ${INFILES_R2}
else
    echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi
