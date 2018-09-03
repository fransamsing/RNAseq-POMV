#!/bin/bash

# Get the md5sum from Downloaded files and check against 
# md5sum provided by Macrogen in the Raw Data Report in Docs/


filenames="$(cut -d , -f 1 ../METADATA.csv | grep -v sample_id)"

for f in $filenames
	do md5sum /OSM/CBR/AF_POMV/work/POMV_RNA_seq/Data/${f} > ../Outputs/md5sum.txt
done
