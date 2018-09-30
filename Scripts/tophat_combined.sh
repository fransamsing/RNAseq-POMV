#!/bin/bash

#***************************************************************#
#              Tophat Alignment with Salmon Genome              #
#                STEP 2 IN BOWTIE-TOPHAT PIPELINE               #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=TOPHAT
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Logs/tophat_%A.out
#SBATCH --array=0-14

module load jemalloc
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE
module load tophat/2.1.1
module load bowtie/2.2.9

# Working Directories
INPDIR=/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Data
#REFDIR=/flush3/sam079/RNAseq-POMV/GenomeIndex/Bowtie_SalmonPOMV
REFDIR=/flush3/sam079/RNAseq-POMV/GenomeIndex/Bowtie_SalmonISAV
#OUTDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment/AlignSalmonPOMVCombined
OUTDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment/AlignSalmonISAVCombined

SAMPLES=( $(cut -d , -f 1 ../STARInputList.csv) );
INFILES_R1_LIST=( $(cut -d , -f 2 ../STARInputList.csv) );
INFILES_R2_LIST=( $(cut -d , -f 3 ../STARInputList.csv) );

if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
then
    i=$SLURM_ARRAY_TASK_ID
    INFILES_R1=${INPDIR}/${INFILES_R1_LIST[$i]}
    INFILES_R2=${INPDIR}/${INFILES_R2_LIST[$i]}
    tophat -p 8 --read-mismatches 10 --read-edit-dist 10 --max-multihits 10 -o ${OUTDIR}/${SAMPLES[$i]} ${REFDIR}/SalmonISAV ${INFILES_R1} ${INFILES_R2}
else
    echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi



## PARAMETERS TWEAKING
## --read-mistmatches default for TOPHAT = 2, set to 10 for better results
## --max-multihits default for TOPHAT and STAR = 20, set to 10 for better results

