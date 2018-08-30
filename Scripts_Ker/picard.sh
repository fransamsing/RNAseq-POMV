#!/bin/bash

#***************************************************************#
#                            picard.sh                          #
#                  written by Kerensa McElroy                   #
#                          May 17, 2018                         #
#                                                               #
#            analyse alignment quality using picard             #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=picard
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=10GB
#SBATCH --output=logs/slurm/mdup_%A_%a.out


#------------------------project variables----------------------#
IN_DIR=${BIG}/analysis/bwa
OUT_DIR=${BIG}/analysis/bwa/picard

#---------------------------------------------------------------#

module add picard/2.9.2 
module add samtools/0.1.19

MeanQualityByCycle --version logs/${TODAY}_main.log

IN_LIST=( $(tail -n +2 ${METADATA} | cut -f 1 | sed "s/_R[12].*//" | sort -u) );

mkdir -p ${OUT_DIR}

if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
    then
    MeanQualityByCycle \
        VALIDATION_STRINGENCY=$VALIDATION \
        R=$REF \
        INPUT=$INDIR/${STEM}.bam \
        OUTPUT=$OUTDIR/${STEM}_MQBC.txt \
        CHART_OUTPUT=$OUTDIR/${STEM}_MQBC.pdf \
        TMP_DIR=/OSM/CBR/NRCA_FINCHGENOM/temp/picard.temp \
        MAX_RECORDS_IN_RAM=$MAX_REC 2> ${OUTDIR}/${STEM}_MQBC.log
    wait

    java -jar -Xmx15g -Djava.io.tmpdir=/OSM/CBR/NRCA_FINCHGENOM/temp /apps/picard/1.138/picard.jar CollectAlignmentSummaryMetrics \
        VALIDATION_STRINGENCY=$VALIDATION \
        LEVEL=$METRIC_LEVEL \
        R=$REF \
        INPUT=$INDIR/${STEM}.bam \
        OUTPUT=$OUTDIR/${STEM}_CASM.txt \
        TMP_DIR=$TMP/picard.temp \
        MAX_RECORDS_IN_RAM=$MAX_REC 2> ${OUTDIR}/${STEM}_CASM.log
    wait

    java -jar -Xmx15g -Djava.io.tmpdir=/OSM/CBR/NRCA_FINCHGENOM/temp /apps/picard/1.138/picard.jar QualityScoreDistribution \
        VALIDATION_STRINGENCY=$VALIDATION \
	STEM=${IN_LIST["$SLURM_ARRAY_TASK_ID"]}
	        VALIDATION_STRINGENCY=$VALIDATION \
        R=$REF \
        INPUT=$INDIR/${STEM}.bam \
        OUTPUT=$OUTDIR/${STEM}_QSD.txt \
        CHART=$OUTDIR/${STEM}_QSD.pdf \
        TMP_DIR=$TMP/picard.temp \
        MAX_RECORDS_IN_RAM=$MAX_REC
    wait

    java -jar -Xmx15g -Djava.io.tmpdir=/OSM/CBR/NRCA_FINCHGENOM/temp /apps/picard/1.138/picard.jar CollectInsertSizeMetrics \
        VALIDATION_STRINGENCY=$VALIDATION \
        R=$REF \
        INPUT=$INDIR/${STEM}.bam \
        OUTPUT=$OUTDIR/${STEM}_CISM.txt \
        HISTOGRAM_FILE=$OUTDIR/${STEM}_CISM.pdf \
        TMP_DIR=$TMP/picard.temp \
        MAX_RECORDS_IN_RAM=$MAX_REC
    wait

    java -jar -Xmx15g -Djava.io.tmpdir=/OSM/CBR/NRCA_FINCHGENOM/temp /apps/picard/1.138/picard.jar CollectWgsMetrics \
        VALIDATION_STRINGENCY=$VALIDATION \
        R=$REF \
        INPUT=$INDIR/${STEM}.bam \
        OUTPUT=$OUTDIR/${STEM}_WGS.txt \
        TMP_DIR=$TMP/picard.temp \
        MAX_RECORDS_IN_RAM=$MAX_REC 2> ${OUTDIR}/${STEM}_WGS.log
    wait

    java -jar -Xmx15g -Djava.io.tmpdir=/OSM/CBR/NRCA_FINCHGENOM/temp /apps/picard/1.138/picard.jar CollectGcBiasMetrics \
        VALIDATION_STRINGENCY=$VALIDATION \
        R=$REF \
        INPUT=$INDIR/${STEM}.bam \
        OUTPUT=$OUTDIR/${STEM}_CGcBM.txt \
        CHART=$OUTDIR/${STEM}_CGcBM.pdf \
        SUMMARY_OUTPUT=$OUTDIR/${STEM}_CGcBM_summary.txt \
        WINDOW_SIZE=$WIN \
        TMP_DIR=$TMP/picard.temp \
        MAX_RECORDS_IN_RAM=$MAX_REC
    wait

    java -jar -Xmx15g -Djava.io.tmpdir=/OSM/CBR/NRCA_FINCHGENOM/temp /apps/picard/1.138/picard.jar EstimateLibraryComplexity \
        VALIDATION_STRINGENCY=$VALIDATION \
        R=$REF \
        INPUT=$INDIR/${STEM}.bam \
        OUTPUT=$OUTDIR/${STEM}_ELC.txt \
        TMP_DIR=$TMP/picard.temp \
        MAX_RECORDS_IN_RAM=$MAX_REC

	
    else
        echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi

mv ${BIG}/logs/slurm/mdup_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out $BIG/logs/${TODAY}_mdup_slurm/


