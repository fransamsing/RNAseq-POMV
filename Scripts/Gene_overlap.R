## Gene gene overlap ## 

setwd("/Volumes/HOME_INTEL/RNAseq-POMV")

source("https://bioconductor.org/biocLite.R")
biocLite("GOplot")
library(GOplot)

POMV6 <- read.csv('Results/ControlvsPOMV6_ALL.csv')
POMV24 <- read.csv('Results/ControlvsPOMV24_ALL.csv')
ISAV6 <- read.csv('Results/ControlvsISAV6_ALL.csv')
ISAV24 <- read.csv('Results/ControlvsISAV24_ALL.csv')

POMV6_SigGeneList <- POMV6[abs(POMV6$logFC) >= 2,]
POMV6_SigGeneList <- POMV6_SigGeneList[,1:2]
write.table(POMV6_SigGeneList, 'Results/POMV6_SigGeneList.csv', sep = ',', row.names = FALSE)

POMV24_SigGeneList <- POMV24[abs(POMV24$logFC) >= 5,]
POMV24_SigGeneList <- POMV24_SigGeneList[,1:2]

ISAV6_SigGeneList <- ISAV6[abs(ISAV6$logFC) >= 2,]
ISAV6_SigGeneList <- ISAV6_SigGeneList[,1:2]
write.table(ISAV6_SigGeneList, 'Results/ISAV6_SigGeneList.csv', sep = ',', row.names = FALSE)

ISAV24_SigGeneList <- ISAV24[abs(ISAV24$logFC) >= 5,]
ISAV24_SigGeneList <- ISAV24_SigGeneList[,1:2]

venn <- GOVenn(POMV6_SigGeneList, ISAV6_SigGeneList, label = c('POMV6' , 'ISAV6'))
venn

intersection_6_UP <- data.frame(intersect(POMV6_SigGeneList$Row.names[POMV6_SigGeneList$logFC > 0], ISAV6_SigGeneList$Row.names[ISAV6_SigGeneList$logFC > 0]))
colnames(intersection_6_UP) <- "SYMBOL"
intersection_6_DOWN <- data.frame(intersect(POMV6_SigGeneList$Row.names[POMV6_SigGeneList$logFC < 0], ISAV6_SigGeneList$Row.names[ISAV6_SigGeneList$logFC < 0]))
colnames(intersection_6_DOWN) <- "SYMBOL"

unique_POMV6_UP <- POMV6_SigGeneList$Row.names[POMV6_SigGeneList$logFC>0][!(POMV6_SigGeneList$Row.names[POMV6_SigGeneList$logFC>0] %in% ISAV6_SigGeneList$Row.name[ISAV6_SigGeneList$logFC>0])]
unique_POMV6_UP <- data.frame(as.character(unique_POMV6_UP))
colnames(unique_POMV6_UP) <- "SYMBOL"

unique_POMV6_DOWN <- POMV6_SigGeneList$Row.names[POMV6_SigGeneList$logFC<0][!(POMV6_SigGeneList$Row.names[POMV6_SigGeneList$logFC<0] %in% ISAV6_SigGeneList$Row.name[ISAV6_SigGeneList$logFC<0])]
unique_POMV6_DOWN <- data.frame(as.character(unique_POMV6_DOWN))
colnames(unique_POMV6_DOWN) <- "SYMBOL"

unique_ISAV6_UP <- ISAV6_SigGeneList$Row.names[ISAV6_SigGeneList$logFC>0][!(ISAV6_SigGeneList$Row.names[ISAV6_SigGeneList$logFC>0] %in% POMV6_SigGeneList$Row.name[POMV6_SigGeneList$logFC>0])]
unique_ISAV6_UP <- data.frame(as.character(unique_ISAV6_UP))
colnames(unique_ISAV6_UP) <- "SYMBOL"

unique_ISAV6_DOWN <- ISAV6_SigGeneList$Row.names[ISAV6_SigGeneList$logFC<0][!(ISAV6_SigGeneList$Row.names[ISAV6_SigGeneList$logFC<0] %in% POMV6_SigGeneList$Row.name[POMV6_SigGeneList$logFC<0])]
unique_ISAV6_DOWN <- data.frame(as.character(unique_ISAV6_DOWN))
colnames(unique_ISAV6_DOWN) <- "SYMBOL"

## Gene Ontology Analysis ##
## get an OrgDb for Atlantic salmon
library(Category)
library(AnnotationHub)
hub <- AnnotationHub()
query(hub, c("salmo salar","orgdb"))
salmodb <- hub[["AH61820"]]
DatPkgFactory(salmodb)
columns(salmodb)

intersection_6_UP_names <- as.character(intersection_6_UP$SYMBOL)
intersection_6_UP_ann <-  AnnotationDbi::select(salmodb, keys = intersection_6_UP_names, columns=c("ENTREZID","GENENAME"), keytype = "SYMBOL")
write.table(intersection_6_UP_ann, 'Results/intersection_6_UP_ann.csv', sep = ',', row.names = FALSE)

intersection_6_DOWN_names <- as.character(intersection_6_DOWN$SYMBOL)
intersection_6_DOWN_ann <-  AnnotationDbi::select(salmodb, keys = intersection_6_DOWN_names, columns=c("ENTREZID","GENENAME"), keytype = "SYMBOL")
write.table(intersection_6_DOWN_ann, 'Results/intersection_6_DOWN_ann.csv', sep = ',', row.names = FALSE)

unique_POMV6_UP_names <- as.character(unique_POMV6_UP$SYMBOL)
unique_POMV6_UP_ann <-  AnnotationDbi::select(salmodb, keys = unique_POMV6_UP_names, columns=c("ENTREZID","GENENAME"), keytype = "SYMBOL")
write.table(unique_POMV6_UP_ann, 'Results/unique_POMV6_UP_ann.csv', sep = ',', row.names = FALSE)

unique_POMV6_DOWN_names <- as.character(unique_POMV6_DOWN$SYMBOL)
unique_POMV6_DOWN_ann <-  AnnotationDbi::select(salmodb, keys = unique_POMV6_DOWN_names, columns=c("ENTREZID","GENENAME"), keytype = "SYMBOL")
write.table(unique_POMV6_DOWN_ann, 'Results/unique_POMV6_DOWN_ann.csv', sep = ',', row.names = FALSE)

unique_ISAV6_UP_names <- as.character(unique_ISAV6_UP$SYMBOL)
unique_ISAV6_UP_ann <-  AnnotationDbi::select(salmodb, keys = unique_ISAV6_UP_names, columns=c("ENTREZID","GENENAME"), keytype = "SYMBOL")
write.table(unique_ISAV6_UP_ann, 'Results/unique_ISAV6_UP_ann.csv', sep = ',', row.names = FALSE)

unique_ISAV6_DOWN_names <- as.character(unique_ISAV6_DOWN$SYMBOL)
unique_ISAV6_DOWN_ann <-  AnnotationDbi::select(salmodb, keys = unique_ISAV6_DOWN_names, columns=c("ENTREZID","GENENAME"), keytype = "SYMBOL")
write.table(unique_ISAV6_DOWN_ann, 'Results/unique_ISAV6_DOWN_ann.csv', sep = ',', row.names = FALSE)
