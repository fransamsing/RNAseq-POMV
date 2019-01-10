# Plot of gene coverage
setwd("/Volumes/sam079/RNAseq-POMV")
library(tidyverse)

#R1 <- read.csv("Data/POMV_coverage/POMV6HPIR1_coverage.txt", sep = '\t', header = FALSE)
#R2 <- read.csv("Data/POMV_coverage/POMV6HPIR2_coverage.txt", sep = '\t', header = FALSE)
#R3 <- read.csv("Data/POMV_coverage/POMV6HPIR3_coverage.txt", sep = '\t', header = FALSE)

#R1 <- read.csv("Data/POMV_coverage/POMV24HPIR1_coverage.txt", sep = '\t', header = FALSE)
#R2 <- read.csv("Data/POMV_coverage/POMV24HPIR2_coverage.txt", sep = '\t', header = FALSE)
#R3 <- read.csv("Data/POMV_coverage/POMV24HPIR3_coverage.txt", sep = '\t', header = FALSE)

#R1 <- read.csv("Data/ISAV_coverage/ISAV6HPIR1_coverage.txt", sep = '\t', header = FALSE)
#R2 <- read.csv("Data/ISAV_coverage/ISAV6HPIR2_coverage.txt", sep = '\t', header = FALSE)
#R3 <- read.csv("Data/ISAV_coverage/ISAV6HPIR2_coverage.txt", sep = '\t', header = FALSE)

R1 <- read.csv("Data/ISAV_coverage/ISAV24HPIR1_coverage.txt", sep = '\t', header = FALSE)
R2 <- read.csv("Data/ISAV_coverage/ISAV24HPIR2_coverage.txt", sep = '\t', header = FALSE)
R3 <- read.csv("Data/ISAV_coverage/ISAV24HPIR2_coverage.txt", sep = '\t', header = FALSE)

newdf <- transpose(list(R1, R2, R3))
segments <- data.frame(newdf$V1)
colnames(segments) <- c("R1", "R2", "R3")
indexes <- data.frame(newdf$V2)
colnames(indexes) <- c("R1", "R2", "R3")
cov <- data.frame(newdf$V3)
colnames(cov) <- c("R1", "R2", "R3")
cov <- cov %>% 
        rowwise() %>% 
        mutate(mean_cov=mean(c(R1, R2, R3)))

cov_ann <- data.frame(cbind(as.character(segments$R1), indexes$R1, cov$mean_cov))
colnames(cov_ann) <- c("segment", "pos", "cov")

## FOR POMV
#cov_ann2 <- cov_ann %>% extract(col = segment, into = "segment", regex = regex("14-01514_POMV_([:alnum:]+)")) 

##FOR ISAV
cov_ann2 <- cov_ann
levels(cov_ann2$segment) <- c("Seg1:PB2", "Seg2:PB1", "Seg3:NP", "Seg4:PA", "Seg5:F", "Seg6:HE", "Seg7:NS1,NEP", "Seg8:M1,M2")

## Fix positions
cov_ann2$pos <- as.numeric(as.character(cov_ann2$pos))
cov_ann2$cov <- as.numeric(as.character(cov_ann2$cov))
str(cov_ann2)

ggplot(cov_ann2, aes(x = pos, y = cov)) + 
  geom_area(aes(fill=segment)) + 
  facet_wrap(~ segment, ncol = 2) +
  labs(x = "Position (bp)", 
       y = "Depth", 
       fill = "Segments")

  




