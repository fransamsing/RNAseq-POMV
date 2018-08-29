#!/bin/bash

#		runs fastqc on gzipped raw data files        
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=fastqc
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Analyses/logs/slurm/%A-%a.out
#SBATCH --array=1-30

#---------------------------------------------------------------#

module load fastqc

cd /OSM/CBR/AF_POMV/work/POMV_RNA_seq/Data/

x=$(sed -n ${SLURM_ARRAY_TASK_ID}p dir_list.txt)
echo $x

fastqc -o ../Analysis/FastQC/ $x


