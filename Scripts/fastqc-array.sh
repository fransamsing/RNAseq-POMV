#!/bin/bash

#***************************************************************#
#		runs fastqc on gzipped raw data files           #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=fastqc
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=slurm_%j.out
#SBATCH --array=1-30

#---------------------------------------------------------------#

module load fastqc

cut -d , -f 1 ../METADATA.csv | grep -v sample_id > /OSM/CBR/AF_POMV/work/POMV_RNA_seq/Data/file_list.txt

cd /OSM/CBR/AF_POMV/work/POMV_RNA_seq/Data/

x=$(sed -n ${SLURM_ARRAY_TASK_ID}p file_list.txt)
echo $x

fastqc -o ../Analysis/FastQC/ $x


