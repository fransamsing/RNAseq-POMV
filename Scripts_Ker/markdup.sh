#!/bin/bash

#***************************************************************#
#                            markdup.sh                         #
#                  written by Kerensa McElroy                   #
#                          May 6, 2018                          #
#                                                               #
#                      mark duplicates using picard             #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=markdup
#SBATCH --time=01:00:00
#SBATCH --ntasks=20
#SBATCH --mem=10GB
#SBATCH --output=logs/slurm/mdup_%A_%a.out


#------------------------project variables----------------------#
IN_DIR=${BIG}/analysis/bwa/orig
OUT_DIR=${BIG}/analysis/bwa/dedup

#---------------------------------------------------------------#

module add picard/2.9.2 
module add samtools

MarkDuplicates --version logs/${TODAY}_main.log
samtools --version logs/${TODAY}_main.log

# mark duplicate reads in sorted bam and merge by sample
IN_LIST=( $(tail -n +2 ${METADATA} | cut -f 2 | sort -u) );

mkdir -P ${OUT_DIR}

if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
    then
	    SAMPLE=${IN_LIST["$SLURM_ARRAY_TASK_ID"]}
            files=`grep "${SAMPLE}" ${METDATA} | cut -f1 | sed "s/_R[12].*/_fixmate_sort.bam/" | sort -u | sed "s_^_I=\.\./orig/_"`
            cd ${OUT_DIR}

            echo '#!/bin/bash' > ${SAMPLE}.picard
            echo "MarkDuplicates \\" >> ${SAMPLE}.picard
            echo $files | xargs -n 1 sh -c 'echo "\t${0} \\" >> ${SAMPLE}.picard'
            echo "  O=${SAMPLE}_fixmate_sort_MDUP.bam \\" >> ${SAMPLE}.picard
            echo "  REMOVE_DUPLICATES=false \\" >> ${SAMPLE}.picard
            echo "  METRICS_FILE=${SAMPLE}_MDUP.txt \\" >> ${SAMPLE}.picard
            echo "  MAX_FILE_HANDLES=1000 \\" >> ${SAMPLE}.picard
            echo "  VALIDATION_STRINGENCY=SILENT 2> ${SAMPLE}_mdup.log "    >> ${SAMPLE}.picard

            chmod a+x ${SAMPLE}.picard
            ${SAMPLE}.picard
            samtools index ${SAMPLE}_fixmate_sort_MDUP.bam 
    else
        echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi

mv ${BIG}/logs/slurm/mdup_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out $BIG/logs/${TODAY}_mdup_slurm/


