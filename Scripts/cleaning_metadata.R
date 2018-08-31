## Open Metadata and clean ## 
# tidyverse 1.2.1
library(tidyverse)


raw_data <- read.csv('Raw_data_file_from_provider.txt', sep = " ", header = TRUE)
raw_data <- raw_data[,1:3]
colnames(raw_data) <- c('sample_id', 'file_size', 'md5sum')
raw_data <- as.tibble(raw_data)

negatives <- raw_data[1:6,1:3]
viruses <- raw_data[7:30,1:3]


viruses_new <- viruses %>% extract(col = sample_id, into = c("treatment", "hpi", "replicate", "read"), 
                             regex = regex("([A-Z]{4})(\\d+).+(\\d)_(\\d)")) 

negatives_new <- negatives %>% extract(col = sample_id, into = c("treatment", "hpi", "replicate", "read"), 
                                 regex = regex("(\\D+)control(R)(\\d)_(\\d)")) %>%
                                mutate_if(is.character, str_replace_all, pattern = "R", replacement = "0")


samples_metadata <- rbind(negatives_new, viruses_new)
samples_metadata

samples_metadata_final <- cbind(raw_data[,1], samples_metadata[,5:6], samples_metadata[,1:4])
samples_metadata_final
str(samples_metadata_final)


write.table(samples_metadata_final, "METADATA.csv", row.names = F, sep=",", quote=FALSE)



