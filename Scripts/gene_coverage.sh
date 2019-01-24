#!/bin/bash

#***************************************************************#
#            Using bedtools to get gene coveraged  		#
#								#
# Gene coverage is scaled by number of reads that aligned 	#
# to each viral genome using samtools flagstats in script 	#
# get_viral_alignment_summary.sh 				#
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=gene_coverage
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Logs/gene_coverage_%A.out
#SBATCH --array=0-14

#-----------------------------------------------------------------#

module load jemalloc
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE
module load samtools/0.1.19
module load bedtools/2.26.0

# Working Directories
#INPDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment/AlignPOMVGenome 
#INPDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment/AlignISAVGenome
INPDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment/AlignPOMVNewGenome 

OUTDIR=/home/sam079/RNAseq-POMV/Data
SAMPLES=( $(cut -d , -f 1 ../STARInputList.csv) );

SCALE=( $(cut -d , -f 4 ../Results/summary_virus_alignment_POMV.csv) )
#SCALE=( $(cut -d , -f 4 ../Results/summary_virus_alignment_ISAV.csv) )

if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
then
    i=$SLURM_ARRAY_TASK_ID
    #genomeCoverageBed -i ${INPDIR}/${SAMPLES[$i]}.sorted.bed -g ../POMV_sequences.txt -d > ${OUTDIR}/POMV_coverage/${SAMPLES[$i]}_coverage.txt -scale ${SCALE[$i]}
	#genomeCoverageBed -i ${INPDIR}/${SAMPLES[$i]}.sorted.bed -g ../ISAV_sequences.txt -d > ${OUTDIR}/ISAV_coverage/${SAMPLES[$i]}_coverage.txt -scale ${SCALE[$i]}
	genomeCoverageBed -i ${INPDIR}/${SAMPLES[$i]}.sorted.bed -g ../POMV_new_sequences.txt -d > ${OUTDIR}/POMV_new_coverage/${SAMPLES[$i]}_coverage.txt -scale ${SCALE[$i]}
else
    echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi

## After running copy files to ~/RNAseq-POMV/Data



