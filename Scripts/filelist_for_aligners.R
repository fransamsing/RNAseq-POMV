# Get METADATA with Read 1 and Read 2 in separate columns for STAR analysis
library(tidyverse)

metadata <- read_csv("METADATA.csv")
metadata

read1 <- metadata %>% filter(read == 1)
read2 <- metadata %>% filter(read == 2)

read1 <- read1[,1]
read2 <- read2[,1]

reads_combined <- cbind(read1, read2)
colnames(reads_combined) <- c('read1', 'read2')

STAR_input_files <- reads_combined %>% separate(read1, c('sample_ID', 'file_format'), sep = '_', remove = FALSE)
STAR_input_files <- cbind(STAR_input_files$sample_ID, STAR_input_files$read1, STAR_input_files$read2)
STAR_input_files 

write.table(STAR_input_files, "STARInputList.csv", row.names = F, col.names = F, sep=",", quote=FALSE)

                             