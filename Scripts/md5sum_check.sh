#!/bin/bash

# Cut the md5sum column from METADATA_final.csv, sort and save as an intermediate file
cut -d , -f3 ../METADATA.csv | grep -v md5sum| sort > ../Results/md5sum_from_provider.csv

# Run md5sum command on data files in Bowen, 
# cut and sort the column with md5sum data and save as an intermediate file
md5sum /OSM/CBR/AF_POMV/work/POMV_RNA_seq/Data/*.gz | cut -d ' ' -f1 | sort > ../Results/md5sum_from_data.csv

# Check the difference in two files

diff ../Results/md5sum_from_provider.csv ../Results/md5sum_from_data.csv

