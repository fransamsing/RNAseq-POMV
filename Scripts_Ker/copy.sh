#!/bin/bash

#***************************************************************#
#                            copy.sh                            #
#                  written by Kerensa McElroy                   #
#                         Feburary 1, 2018                      #
#                                                               #
#               batch file copy including integrity check       #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=copy
#SBATCH --partition=io
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=1GB
#SBATCH --output=logs/slurm/copy_%A_%a.out

#-------------------------commands-------------------------#

IN_FILE_LIST=(${IN_DIR}/*${IN_EXT})

if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
    then
        IN_FILE=${IN_FILE_LIST["$SLURM_ARRAY_TASK_ID"]}
        BASE=$(basename ${IN_FILE})
        md5sum ${IN_FILE} > $OUT_DIR/${BASE}.md5
        cp -n --no-preserve=ownership ${IN_FILE} $OUT_DIR/${BASE}
        md5sum -c $OUT_DIR/${BASE}.md5 >> $BIG/logs/${TODAY}_copy.log
    else
        echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi

mkdir -p logs/${TODAY}_copy_slurm
mv logs/slurm/copy_* logs/${TODAY}_copy_slurm/
