#!/bin/bash

#***************************************************************#
#                         bwa_index.sh                          #
#                  written by Kerensa McElroy                   #
#                         May 3, 2018                           #
#                                                               #
#                   generate indexed reference                  #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=bwa_index
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --mem=5GB
#SBATCH --output=logs/slurm/bwa_index.out


#------------------------project variables----------------------#
IN_DIR=${BIG}/data

#---------------------------------------------------------------#

module add bwa

bwa | head -5 >> ${BIG}/logs/${TODAY}_main.log
echo '' >> ${BIG}/logs/${TODAY}_main.log

bwa index -a bwtsw -p ${IN_DIR}/${REF%.*}  ${IN_DIR}/${REF} \
    2>> ${BIG}/logs/${TODAY}_main.log


mv ${BIG}/logs/slurm/bwa_index.out $BIG/logs/${TODAY}_bwa_index_slurm/


