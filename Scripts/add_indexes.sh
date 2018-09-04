#!/bin/bash

# To create a new .csv file with indexes
for f in $filenames; do zcat /OSM/CBR/AF_POMV/work/POMV_RNA_seq/Data/${f} | head -n 1 ; done |cut -f10 -d ':' | tr '+' ',' > index_columns.csv

# To add a header to above column which will create a new column
echo -e "index1,tindex2" | cat - index_columns.csv > index_columns1.csv

# To combine the index_columns1.csv with METADATA.csv

paste -d, METADATA.csv index_columns1.csv > METADATA1.csv
