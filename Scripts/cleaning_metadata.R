## Open Metadata and clean ## 

# tidyverse 1.2.1
library(tidyverse)

samples <- read.table('metadata.txt', sep = " ")
samples <- as.tibble(samples)
samples
colnames(samples) <- c('sample_id','total_read_bases_(bp)','total_reads', 'GC','AT','Q20','Q30')

## DATA WITHOUT SAMPLE IDs ##
samples_data <- samples[,2:7]
samples_data

## SAMPLE IDs without Neg controls
samples_id_virus <- samples[4:15,1]
samples_id_virus

samples_id_virus<- samples_id_virus %>% extract(col = sample_id, into = c("treatment", "hpi", "replicate"), 
                             regex = regex("(\\D+)(\\d+).+(\\d)"))

# Neg controls IDs 
samples_id_neg <- samples[1:3,1]
samples_id_neg <- samples_id_neg %>% extract(col = sample_id, into = c("treatment", "hpi", "replicate"), 
                          regex = regex("(\\D+)(.+)(\\d)")) %>% 
                  mutate_if(is.character, str_replace_all, pattern = "R", replacement = "0")

samples_id_data <- rbind(samples_id_neg, samples_id_virus)
samples_id_data

samples_unique_id <- samples[,1]
samples_metadata <- cbind(samples_unique_id, samples_id_data, samples_data)   
samples_metadata <- as.tibble(samples_metadata)
samples_metadata

write.csv(samples_metadata, "METADATA", row.names = FALSE)
