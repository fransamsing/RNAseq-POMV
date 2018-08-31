#!/bin/bash


#***************************************************************#
#                        get_metadata.sh                        #
#                     written by Fran Samsing                   #
#                        August 31, 2018                        #
#                                                               #
#           get metadata from sequencing files and store        #
#***************************************************************#



filenames="$(cut -d , -f 1 METADATA.csv | grep -v sample_id)"

for f in $filenames;do echo $f; zcat /OSM/CBR/AF_POMV/work/POMV_RNA_seq/Data/*.fastq.gz | head -n 1 ; done
