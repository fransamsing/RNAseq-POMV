#!/bin/bash

#***************************************************************#
#              Merge Transcripts with Cuffmerge                 #
#           STEP 4 IN BOWTIE-TOPHAT-CUFFLINKS PIPELINE          #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=CUFFMERGE
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=15GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Logs/cuffmerge_%A.out


module load cufflinks

# Working Directories

REFDIR=/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Genomes
INPDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment/AlignSalmonPOMVCombined
OUTDIR=/flush3/sam079/RNAseq-POMV/Processed/Assembly/SalmonPOMV/Cufflinks

#ls `find $OUTDIR/NEG* -type f` | grep transcript > $OUTDIR/NEG.txt
#ls `find $OUTDIR/*6* -type f` | grep transcript | sort -r > $OUTDIR/6HPI.txt
#ls `find $OUTDIR/*24* -type f` | grep transcript | sort -r > $OUTDIR/24HPI.txt

#cat $OUTDIR/NEG.txt $OUTDIR/6HPI.txt $OUTDIR/24HPI.txt > $OUTDIR/assemblies.txt

#cuffmerge -s $REFDIR/SalmonPOMV.fa -p 8 -o $OUTDIR/cuffmerge $OUTDIR/assemblies.txt

cuffdiff -o $OUTDIR/cuffdiff -b $REFDIR/SalmonPOMV.fa -p 8 -L NEG,POMV6HPI,ISAV6HPI,POMV24HPI,ISAV24HPI \
-u $OUTDIR/cuffmerge/merged.gtf \
$INPDIR/NEGcontrolR1/accepted_hits.bam,$INPDIR/NEGcontrolR2/accepted_hits.bam,$INPDIR/NEGcontrolR3/accepted_hits.bam \
$INPDIR/POMV6HPIR1/accepted_hits.bam,$INPDIR/POMV6HPIR2/accepted_hits.bam,$INPDIR/POMV6HPIR3/accepted_hits.bam \
$INPDIR/ISAV6HPIR1/accepted_hits.bam,$INPDIR/ISAV6HPIR2/accepted_hits.bam,$INPDIR/ISAV6HPIR3/accepted_hits.bam \
$INPDIR/POMV24HPIR1/accepted_hits.bam,$INPDIR/POMV24HPIR2/accepted_hits.bam,$INPDIR/POMV24HPIR3/accepted_hits.bam \
$INPDIR/ISAV24HPIR1/accepted_hits.bam,$INPDIR/ISAV24HPIR2/accepted_hits.bam,$INPDIR/ISAV24HPIR3/accepted_hits.bam
 

