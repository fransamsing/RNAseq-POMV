#!/bin/bash

#***************************************************************#
#                Align Raw Reads Against Salmon Genome          #
#                      STEP 2 IN STAR PIPELINE                  #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=STAR_align
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Logs/star_%A.out
#SBATCH --array=1-15


module load jemalloc
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE
module load star

# Working Dictories 
INPDIR=/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Data
ANODIR=/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Genomes/Salmo_salar
REFDIR=/flush3/sam079/RNAseq-POMV/GenomeIndex
OUTDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment

SAMPLES=( $(cut -d , -f 1 ../STARInputList.csv) );
INFILES_R1=( $(cut -d , -f 2 ../STARInputList.csv) );
INFILES_R2=( $(cut -d , -f 3 ../STARInputList.csv) );


if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
then
    i=$SLURM_ARRAY_TASK_ID
    STAR \
    --runThreadN 8 \
    --genomeDir ${REFDIR}/ \
    --readFilesIn ${INPDIR}/${INFILES_R1[$i]} ${INPDIR}/${INFILES_R2[$i]} \
    --readFilesCommand zcat \
    --outFileNamePrefix ${OUTDIR}/${SAMPLES[$i]} \
    --outFilterMismatchNmax 10 \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM GeneCounts \
    --outReadsUnmapped Fastx \
    --outFilterMultimapNmax 10 \
    --outSAMattrIHstart 0 \
    --outSAMmapqUnique 255 \
    --outSAMmultNmax -1 \
    --chimSegmentMin 40
else
    echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi
