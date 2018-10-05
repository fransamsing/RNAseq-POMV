#!/bin/bash

#***************************************************************#
#                  Estimate Transcript Abundances               #
#                STEP 5 IN HISAT-STRINGTIE PIPELINE             #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=Stringtie_ballgown
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Logs/stringtie_ballgown_%A.out
#SBATCH --array=0-14

module load jemalloc
export OMP_NUM_THREADS=$SLURM_NTASKS_PER_NODE
module load stringtie

# Working Directories
REFDIR=/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Genomes
ANODIR=/flush3/sam079/RNAseq-POMV/Processed/Assembly/SalmonPOMV
INPDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment/AlignSalmonPOMVHisat
OUTDIR=/flush3/sam079/RNAseq-POMV/Processed/Expression/SalmonPOMV/Transcript_abundance

SAMPLES=( $(cut -d , -f 1 ../STARInputList.csv) );


if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
then
    i=$SLURM_ARRAY_TASK_ID
    INFILES=${INPDIR}/${SAMPLES[$i]}
    OUTFILES=${OUTDIR}/${SAMPLES[$i]}
    stringtie -e -B -p 8 -G ${ANODIR}/SalmonPOMV_stringtie_merged.gtf -o ${OUTFILES}/transcripts.gtf -A ${OUTFILES}/gene_abundances.tsv ${INFILES}.bam
else
    echo "Error: Missing array index as SLURM_ARRAY_TASK_ID"
fi

## create filepath lists
ls $OUTDIR > /flush3/sam079/RNAseq-POMV/Processed/Expression/SalmonPOMV/sample_ids.txt
ls -d -1 $OUTDIR/** > /flush3/sam079/RNAseq-POMV/Processed/Expression/SalmonPOMV/filepathlist.txt
paste -d$'\t'  /flush3/sam079/RNAseq-POMV/Processed/Expression/SalmonPOMV/sample_ids.txt /flush3/sam079/RNAseq-POMV/Processed/Expression/SalmonPOMV/filepathlist.txt >> /flush3/sam079/RNAseq-POMV/Processed/Expression/SalmonPOMV/sampleIDs_filepaths.txt

# Create list for the python program 
#ls `find $OUTDIR/NEG* -type f` | grep transcript > $OUTDIR/NEG.txt
#ls `find $OUTDIR/*6* -type f` | grep transcript | sort -r > $OUTDIR/6HPI.txt
#ls `find $OUTDIR/*24* -type f` | grep transcript | sort -r > $OUTDIR/24HPI.txt

#cat $OUTDIR/NEG.txt $OUTDIR/6HPI.txt $OUTDIR/24HPI.txt > $OUTDIR/filepathlist.txt

#paste -d ' '  $OUTDIR/sample_id.txt $OUTDIR/filepathlist.txt >> ../transcript_filepaths.txt
