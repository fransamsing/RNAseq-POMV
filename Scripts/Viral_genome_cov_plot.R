# Plot of gene coverage
setwd("/Volumes/HOME_INTEL/RNAseq-POMV")
library(tidyverse)

POMVR1 <- read.csv("Data/POMV24HPIR1_coverage.txt", sep = '\t', header = FALSE)
POMVR2 <- read.csv("Data/POMV24HPIR2_coverage.txt", sep = '\t', header = FALSE)
POMVR3 <- read.csv("Data/POMV24HPIR3_coverage.txt", sep = '\t', header = FALSE)

#POMVR1 <- read.csv("Data/POMV6HPIR1_coverage.txt", sep = '\t', header = FALSE)
#POMVR2 <- read.csv("Data/POMV6HPIR2_coverage.txt", sep = '\t', header = FALSE)
#POMVR3 <- read.csv("Data/POMV6HPIR3_coverage.txt", sep = '\t', header = FALSE)

newdf <- transpose(list(POMVR1, POMVR2, POMVR3))
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

cov_ann2 <- cov_ann %>% extract(col = segment, into = "segment", regex = regex("14-01514_POMV_([:alnum:]+)")) 
cov_ann2$pos <- as.numeric(as.character(cov_ann2$pos))
cov_ann2$cov <- as.numeric(as.character(cov_ann2$cov))
str(cov_ann2)

ggplot(cov_ann2, aes(x = pos, y = cov)) + 
  geom_area(aes(fill=segment)) + 
  facet_wrap(~ segment, ncol = 2) +
  labs(x = "Position", 
       y = "Depth", 
       fill = "Segments")

  




