#!/bin/bash

#***************************************************************#
#                            bam.sh                             #
#                  written by Kerensa McElroy                   #
#                          May 6, 2018                          #
#                                                               #
#                       sam to bam conversion                   #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=bam
#SBATCH --time=01:00:00
#SBATCH --ntasks=20
#SBATCH --mem=10GB
#SBATCH --output=logs/slurm/bam_%A_%a.out


#------------------------project variables----------------------#
IN_DIR=${BIG}/analysis/bwa

#---------------------------------------------------------------#

module add samtools/1.3.1

samtools --version logs/${TODAY}_main.log

# convert sam to sorted bam format
IN_LIST=( $(tail -n +2 ${METADATA} | cut -f 1 | sed "s/_R[12].*//" | sort -u) );

if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
    then
	    STEM=${IN_LIST["$SLURM_ARRAY_TASK_ID"]}
            cd ${IN_DIR}
            samtools fixmate -O bam ${STEM}*sam ${STEM}_fixmate.bam && rm ${STEM}*sam
            samtools sort -O bam -@ 1 -m 8G -o ${STEM}_fixmate_sort.bam ${STEM}_fixmate.bam && rm ${STEM}_fixmate.bam
            samtools index ${STEM}_fixmate_sort.bam
    else
        echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi

mv ${BIG}/logs/slurm/bam_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out $BIG/logs/${TODAY}_bam_slurm/


