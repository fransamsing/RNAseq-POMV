#!/bin/bash

#***************************************************************#
#                        BOWTIE ALIGN FOR VIRUS                 #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=BOWTIE_ALIGN
#SBATCH --time=00:45:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Logs/bowtie_align_%A.out
#SBATCH --array=0-14

module load jemalloc
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE
module load bowtie/2.2.9

# Working Directories
INPDIR=/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Data
REFDIR=/flush3/sam079/RNAseq-POMV/ViralGenomeIndex/POMVGenomeIndex
OUTDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment/AlignPOMVGenome

SAMPLES=( $(cut -d , -f 1 ../STARInputList.csv) );
INFILES_R1_LIST=( $(cut -d , -f 2 ../STARInputList.csv) );
INFILES_R2_LIST=( $(cut -d , -f 3 ../STARInputList.csv) );

if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
then
    i=$SLURM_ARRAY_TASK_ID
    INFILES_R1=${INPDIR}/${INFILES_R1_LIST[$i]}
    INFILES_R2=${INPDIR}/${INFILES_R2_LIST[$i]}
    bowtie2 -p 8 -x ${REFDIR}/POMV -1 ${INFILES_R1} -2 ${INFILES_R2} -S ${OUTDIR}/${SAMPLES[$i]}.sam 
else
    echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi


