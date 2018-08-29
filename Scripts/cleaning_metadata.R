## Open Metadata and clean ## 

# tidyverse 1.2.1
library(tidyverse)

files <- list.files(path = '/OSM/CBR/AF_POMV/work/POMV_RNA_seq/Data/')
files <- files[1:30]

files <- as.tibble(files)
colnames(files) <- 'sample_id'

negatives <- files[13:18,1]
pomv <- files[19:30,1]
isav <- files[1:12,1]


pomv_new <- pomv %>% extract(col = sample_id, into = c("treatment", "hpi", "replicate", "direction"), 
                             regex = regex("([A-Z]{4})(\\d+).+(\\d)_(\\d)"))

isav_new <- isav %>% extract(col = sample_id, into = c("treatment", "hpi", "replicate", "direction"), 
                             regex = regex("([A-Z]{4})(\\d+).+(\\d)_(\\d)"))

negatives_new <- negatives %>% extract(col = sample_id, into = c("treatment", "hpi", "replicate", "direction"), 
                                 regex = regex("(\\D+)(R)(\\d)_(\\d)")) %>%
                                mutate_if(is.character, str_replace_all, pattern = "R", replacement = "0")


samples_metadata <- rbind(pomv_new, isav_new, negatives_new)
samples_metadata <- cbind(files, samples_metadata)
samples_metadata$hpi <- as.integer(samples_metadata$hpi)
samples_metadata$replicate <- as.integer(samples_metadata$replicate)
samples_metadata$direction <- as.integer(samples_metadata$direction)
str(samples_metadata)

write.csv(samples_metadata, "METADATA.csv", row.names = FALSE)
