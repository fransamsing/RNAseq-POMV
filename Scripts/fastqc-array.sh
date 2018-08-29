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
#SBATCH --output=../Analyses/logs/slurm/*.out

#---------------------------------------------------------------#

module load fastqc

x=$(sed -n ${SLURM_ARRAY_TASK_ID}p dir_list.txt)
echo $x

rsync -a -e "ssh -i $HOME/.ssh/id_rsa_blah" blah@somewhere.nci.org.au:/g/data3/results/fastq/$x BAM_gVCF


for i in /OSM/CBR/AF_POMV/work/POMV_RNA_seq/Data/*.fastq.gz
do
	fastqc -o ../Analyses/FastQC/ $i
done

     
   

