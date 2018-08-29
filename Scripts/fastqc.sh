#!/bin/bash
                            
		                                               
#             runs fastqc on gzipped raw data files             #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=fastqc
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=1GB
#SBATCH --output=logs/slurm/fastqc_%A_%a.out


#------------------------project variables----------------------#
IN_DIR=/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Data
OUT_DIR=/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Results

#---------------------------------------------------------------#

module add fastqc
fastqc -v >> logs/${TODAY}_main.log
mkdir -p $OUT_DIR

IN_FILE_LIST=( $(cut -f 1 ${METADATA}) );

if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
    then
        IN_FILE=${IN_DIR}/${IN_FILE_LIST["$SLURM_ARRAY_TASK_ID"]}
        CMD="fastqc ${IN_FILE} --noextract -f fastq -t 1 -o ${OUT_DIR}"
	echo -e "$CMD"
	${CMD}
    else
        echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi

mv ${BIG}/logs/slurm/fastqc_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out $BIG/logs/${TODAY}_fastqc_slurm/
