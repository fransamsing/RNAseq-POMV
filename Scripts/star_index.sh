 #!/bin/bash

#***************************************************************#
#    		 Generate the Salmon Genome Index               #
#                    STEP 1 IN STAR PIPELINE   	                #
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=STAR_index
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Logs/star_%A.out

module load star

# Working Dictories 
FILE_DIR=/OSM/CBR/AF_POMV/work/POMV_RNA_seq
OUT_DIR=/flush3/sam079/RNAseq-POMV/GenomeIndex


STAR \
--runMode genomeGenerate \
--genomeDir ${OUT_DIR}/Star \
--runThreadN 8 \
--genomeChrBinNbits 14 \
--genomeFastaFiles  ${FILE_DIR}/Genomes/Salmo_salar/GCF_000233375.1_ICSASG_v2_genomic.fna \
--sjdbGTFfile ${FILE_DIR}/Genomes/Salmo_salar/GCF_000233375.1_ICSASG_v2_genomic.gff \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFtagExonParentGene gene \
--sjdbOverhang 100 \
--sjdbGTFfeatureExon exon \
--limitGenomeGenerateRAM=168253832576 \
--outFileNamePrefix ../Logs/star_ 

## --genomeChrBinNbits 14 because 
# log2(GenomeLength/NumberOfReferences) = 13.64 
# Genome Length = 2966890203 (grep -v '^>' GCF_000233375.1_ICSASG_v2_genomic.fna | tr -d '\n' | wc -c)
# NumberOfReferences = 232155 (grep -c '^>' GCF_000233375.1_ICSASG_v2_genomic.fna)


