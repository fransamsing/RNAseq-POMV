#!/bin/bash

#***************************************************************#
#                            trim.sh                            #
#                  written by Kerensa McElroy                   #
#                         March 15, 2018                        #
#                                                               #
#                 trim adapters from raw data                   #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=fastqc
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=5GB
#SBATCH --output=logs/slurm/trim_%A_%a.out


#------------------------project variables----------------------#
IN_DIR=${BIG}/data
OUT_DIR=${BIG}/analysis/trim
ADAPTERS=~/adapters

#---------------------------------------------------------------#

module add trimmomatic

echo "trimmomatic PE version:" >> logs/${TODAY}_main.log
trimmomatic PE -version >> logs/${TODAY}_main.log

mkdir -p $OUT_DIR

IN_LIST=( $(cut -f 1 ${METADATA} | sed "s/_R[12].*//" | sort -u) );

if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
    then
	STEM=${IN_LIST["$SLURM_ARRAY_TASK_ID"]}
        ADAPT_SEQ=`echo ${STEM} | sed "s/.*${ADAPT_LEFT}//;s/${ADAPT_RIGHT}.*//"`
        ADAPT_FILE=`echo ${ADAPTERS}/${ADAPT_TYPE}*${ADAPT_SEQ}*.fa` 
        CMD="trimmomatic PE \
            -threads 1 \
            -phred${PHRED} \
            -trimlog ${OUT_DIR}/${STEM}_log.txt \
            `echo ${IN_DIR}/${STEM}*${EXT}` \
            ${OUT_DIR}/${STEM}_R1_p.fq \
            ${OUT_DIR}/${STEM}_R1_u.fq \
            ${OUT_DIR}/${STEM}_R2_p.fq \
            ${OUT_DIR}/${STEM}_R2_u.fq \
            ILLUMINACLIP:${ADAPT_FILE}:${SEEDMISMATCH}:${PALINCLIP}:${SIMPLECLIP}:${MINADAPTLEN}:${KEEPREADS} \
            SLIDINGWINDOW:${WINDOWSIZE}:${BASEQUAL} \
            LEADING:${BASEQUAL} \
            TRAILING:${BASEQUAL} \
            MINLEN:${MINLENGTH}" 
        echo -e "$CMD" > ${OUT_DIR}/${STEM}.txt
	gzip ${OUT_DIR}/${STEM}*
        ${CMD} >> ${OUT_DIR}/${STEM}.txt 2>&1
    else
        echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi

mv ${BIG}/logs/slurm/trim_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out $BIG/logs/${TODAY}_trim_slurm/


