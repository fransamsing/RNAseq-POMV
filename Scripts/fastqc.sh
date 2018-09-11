#!/bin/bash

#***************************************************************#
#               runs fastqc on gzipped raw data files           #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=fastqc
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Logs/fastqc_%A.out
#SBATCH --array=0-29

#----------------------project variables------------------------# 
IN_DIR=/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Data
OUT_DIR=/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Analysis/FastQC

#---------------------------------------------------------------#

module load fastqc

IN_FILE_LIST=( $(cut -d , -f 1 ../METADATA.csv | grep -v sample_id) );

if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
    then
        i=$SLURM_ARRAY_TASK_ID
        IN_FILE=${IN_DIR}/${IN_FILE_LIST[$i]}
        fastqc ${IN_FILE} --noextract -o ${OUT_DIR}
    else
        echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi

