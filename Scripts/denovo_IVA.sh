#!/bin/bash

#***************************************************************#
#                De novo assembly of POMV using IVA		  	    #                 
#***************************************************************#

#--------------------------sbatch header------------------------#

#SBATCH --job-name=IVA_denovo_assembly
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=francisca.samsingpedrals@csiro.au
#SBATCH --output=../Logs/IVA_%A.out

module load python/3.6.1 
module load kmc 
module load samtools/1.9.0 
module load smalt 
module load mummer

INPDIR=/flush3/sam079/RNAseq-POMV/Processed/Alignment/AlignSalmonGenome/Star
OUTDIR=/flush3/sam079/RNAseq-POMV/Processed/DenovoAssemblyPOMV

~/.local/bin/iva --threads 8 -f ${INPDIR}/POMV24HPIR1Unmapped.out.mate1.fastq -r ${INPDIR}/POMV24HPIR1Unmapped.out.mate2.fastq ${OUTDIR}/IVA


