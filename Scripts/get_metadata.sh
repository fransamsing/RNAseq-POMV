#!/bin/bash


#***************************************************************#
#                        get_metadata.sh                        #
#                     written by Fran Samsing                   #
#                        August 31, 2018                        #
#                                                               #
#           get metadata from sequencing files and store        #
#***************************************************************#


# This script gets the first field of the METADATA.csv file that
# containes the filenames of the RNAseq-samples
# and uses those names to loop through the files stored in the 
# POMV_RNA_seq folder in the OSM storage 

filenames="$(cut -d , -f 1 ../METADATA.csv | grep -v sample_id)"

for f in $filenames;do zcat /OSM/CBR/AF_POMV/work/POMV_RNA_seq/Data/${f} | head -n 1 ; done
