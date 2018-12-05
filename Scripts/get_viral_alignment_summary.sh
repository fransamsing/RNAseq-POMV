#!/bin/bash

module load samtools

FILENAMES=$(cut -d , -f 1 ../STARInputList.csv)
INPDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment/AlignPOMVGenome

for f in $FILENAMES
	do 	echo $f \
		samtools flagstat $INPDIR/${f}.sorted.bam
done > ../Results/summary_virus_alignment.txt
