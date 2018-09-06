# Cut the md5sum column from METADATA_final.csv, sort and save as an intermediate file
cut -d , -f3 METADATA_final.csv | grep -v md5sum| sort > Outputs/METADATA_md5sum.csv

# Run md5sum command on data files in Bowen, 
# cut and sort teh column with md5sum data and save as an intermediate file
md5sum /OSM/CBR/AF_POMV/work/POMV_RNA_seq/Data/*.gz | cut -d ' ' -f1 | sort > Outputs/Data_md5sum.csv

# Check the difference in two files

diff Outputs/METADATA_md5sum.csv Outputs/Data_md5sum.csv 
