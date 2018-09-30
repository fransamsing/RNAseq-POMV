#!/bin/bash
filenames="$(cut -d , -f 1 ../METADATA.csv | grep -v sample_id)"

# To create a new .csv file with indexes
for f in $filenames; do zcat /OSM/CBR/AF_POMV/work/POMV_RNA_seq/Data/${f} | head -n 1 ; done |cut -f10 -d ':' | tr '+' ',' > tmp_indexes.csv

# To add a header to above column which will create a new column
echo -e "index1,index2" | cat - tmp_indexes.csv > ../Results/index_columns.csv

rm tmp_indexes.csv

# To combine the index_columns1.csv with METADATA.csv
metadata=../METADATA.csv

paste -d, $metadata ../Results/index_columns.csv >> ../METADATA_Idx.csv
