#!/bin/bash
module load gmap/2017-11-15 
module load samtools
module load bowtie

REFDIR=/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Genomes/POMV_genome_AAHL 
INPDIR=/flush3/sam079/RNAseq-POMV/Processed/DenovoAssemblyPOMV
READSDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment/AlignSalmonGenome/Star
OUTDIR=/flush3/sam079/RNAseq-POMV/Processed/DenovoAssemblyPOMV

## POMV ORFs
IDXDIR=/flush3/sam079/RNAseq-POMV/ViralGenomeIndex/POMVGenomeIndex

bowtie2 -p 8 -x ${IDXDIR}/POMV -1 ${READSDIR}/POMV24HPIR1Unmapped.out.mate1.fastq \
-2 ${READSDIR}/POMV24HPIR1Unmapped.out.mate2.fastq  \
-S ${OUTDIR}/BOWTIE_unmapped.sam 2> ${OUTDIR}/BOWTIE_unmapped.log 

samtools view -Sb ${OUTDIR}/BOWTIE_unmapped.sam > ${OUTDIR}/BOWTIE_unmapped.bam \

samtools sort ${OUTDIR}/BOWTIE_unmapped.bam -o ${OUTDIR}/BOWTIE_unmapped.sorted.bam \

samtools index ${OUTDIR}/BOWTIE_unmapped.sorted.bam

## POMV New Genome
IDXDIR=/flush3/sam079/RNAseq-POMV/ViralGenomeIndex/POMVNewGenomeIndex

bowtie2 -p 8 -x ${IDXDIR}/POMV -1 ${READSDIR}/POMV24HPIR1Unmapped.out.mate1.fastq \
-2 ${READSDIR}/POMV24HPIR1Unmapped.out.mate2.fastq  \
-S ${OUTDIR}/BOWTIE_new_unmapped.sam 2> ${OUTDIR}/BOWTIE_new_unmapped.log 

samtools view -Sb ${OUTDIR}/BOWTIE_new_unmapped.sam > ${OUTDIR}/BOWTIE_new_unmapped.bam \

samtools sort ${OUTDIR}/BOWTIE_new_unmapped.bam -o ${OUTDIR}/BOWTIE_new_unmapped.sorted.bam \

samtools index ${OUTDIR}/BOWTIE_new_unmapped.sorted.bam

## GMAP INDEX
gmap_build -d POMV_genome -D ${OUTDIR}/GMAP \
-k 13 ${REFDIR}/POMV_14_01514_ORF.fa 

## SPADES
gmap -n 0 -D ${OUTDIR}/GMAP \
-d POMV_genome ${INPDIR}/SPAdes/transcripts.fasta \
-f samse > ${OUTDIR}/GMAP/spades_gmap.sam

samtools view -Sb ${OUTDIR}/GMAP/spades_gmap.sam > ${OUTDIR}/GMAP/spades_gmap.bam \

samtools sort ${OUTDIR}/GMAP/spades_gmap.bam -o ${OUTDIR}/GMAP/spades_gmap.bam \

samtools index ${OUTDIR}/GMAP/spades_gmap.bam

## TRINITY
gmap -n 0 -D ${OUTDIR}/GMAP \
-d POMV_genome ${INPDIR}/Trinity/Trinity.fasta \
-f samse > ${OUTDIR}/GMAP/trinity_gmap.sam

samtools view -Sb ${OUTDIR}/GMAP/trinity_gmap.sam > ${OUTDIR}/GMAP/trinity_gmap.bam \

samtools sort ${OUTDIR}/GMAP/trinity_gmap.bam -o ${OUTDIR}/GMAP/trinity_gmap.bam \

samtools index ${OUTDIR}/GMAP/trinity_gmap.bam

## GMAP NEW INDEX
gmap_build -d POMV_new_genome -D ${OUTDIR}/GMAP \
-k 13 ${REFDIR}/POMV_14_01514_GENOME_format.fa

## SPADES - NEW POMV GENOME
gmap -n 0 -D ${OUTDIR}/GMAP \
-d POMV_new_genome ${INPDIR}/SPAdes/transcripts.fasta \
-f samse > ${OUTDIR}/GMAP/spades_new_gmap.sam

samtools view -Sb ${OUTDIR}/GMAP/spades_new_gmap.sam > ${OUTDIR}/GMAP/spades_new_gmap.bam \

samtools sort ${OUTDIR}/GMAP/spades_new_gmap.bam -o ${OUTDIR}/GMAP/spades_new_gmap.bam \

samtools index ${OUTDIR}/GMAP/spades_new_gmap.bam

## TRINITY - NEW POMV GENOME
gmap -n 0 -D ${OUTDIR}/GMAP \
-d POMV_new_genome ${INPDIR}/Trinity/Trinity.fasta \
-f samse > ${OUTDIR}/GMAP/trinity_new_gmap.sam

samtools view -Sb ${OUTDIR}/GMAP/trinity_new_gmap.sam > ${OUTDIR}/GMAP/trinity_new_gmap.bam \

samtools sort ${OUTDIR}/GMAP/trinity_new_gmap.bam -o ${OUTDIR}/GMAP/trinity_new_gmap.bam \

samtools index ${OUTDIR}/GMAP/trinity_new_gmap.bam

