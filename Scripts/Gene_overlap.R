## Gene gene overlap ## 
setwd("/Volumes/sam079/RNAseq-POMV")

#source("https://bioconductor.org/biocLite.R")
#biocLite("GOplot")
library(GOplot)
library(VennDiagram)
library(tidyverse)

POMV6 <- read.csv('Results/ControlvsPOMV6_ALL.csv')
POMV24 <- read.csv('Results/ControlvsPOMV24_ALL.csv')
ISAV6 <- read.csv('Results/ControlvsISAV6_ALL.csv')
ISAV24 <- read.csv('Results/ControlvsISAV24_ALL.csv')

# Venn Diagrams
POMV6_SigGeneList <- POMV6[abs(POMV6$logFC) >= 2,]
POMV6_SigGeneList <- POMV6_SigGeneList[,c(8,2)]

POMV24_SigGeneList <- POMV24[abs(POMV24$logFC) >= 2,]
POMV24_SigGeneList <- POMV24_SigGeneList[,c(8,2)]
POMV24_SigGeneList <- drop_na(POMV24_SigGeneList)

ISAV6_SigGeneList <- ISAV6[abs(ISAV6$logFC) >= 2,]
ISAV6_SigGeneList <- ISAV6_SigGeneList[,c(8,2)]

ISAV24_SigGeneList <- ISAV24[abs(ISAV24$logFC) >= 2,]
ISAV24_SigGeneList <- ISAV24_SigGeneList[,c(8,2)]

venn6 <- GOVenn(POMV6_SigGeneList, ISAV6_SigGeneList, label = c('POMV6' , 'ISAV6'), plot = FALSE) ## plot = FALSE to get DATA
venn6$plot
unique_POMV6 <- rownames_to_column(data.frame(venn6$table$A_only), var = "ENTREZID")
unique_ISAV6 <- rownames_to_column(data.frame(venn6$table$B_only), var = "ENTREZID")
intersection_6 <- rownames_to_column(data.frame(venn6$table$AB), var = "ENTREZID")

venn24 <- GOVenn(POMV24_SigGeneList, ISAV24_SigGeneList, label = c('POMV24' , 'ISAV24'), plot = FALSE)
venn24$plot
unique_POMV24 <- rownames_to_column(data.frame(venn24$table$A_only), var = "ENTREZID")
unique_ISAV24 <- rownames_to_column(data.frame(venn24$table$B_only), var = "ENTREZID")
intersection_24 <- rownames_to_column(data.frame(venn24$table$AB), var = "ENTREZID")

## DRAW VENN DIAGRAMS 
nPOMV6_UP <- length(unique_POMV6$Trend[unique_POMV6$Trend == "UP"]) + length(intersection_6$Trend[intersection_6$Trend == "UP"])
nISAV6_UP <- length(unique_ISAV6$Trend[unique_ISAV6$Trend == "UP"]) + length(intersection_6$Trend[intersection_6$Trend == "UP"])
nINT6_UP <- length(intersection_6$Trend[intersection_6$Trend == "UP"])

UP6 <- draw.pairwise.venn(nPOMV6_UP, nISAV6_UP, nINT6_UP, category = c("POMV 6 hpi", "ISAV 6 hpi"), 
                   lty = rep("blank", 2), fill = c("red", "orange"), alpha = rep(0.5, 2), cat.pos = c(170, 180), cat.dist = c(0.006, 0.006), scaled = TRUE, rotation.degree = 180, 
                   cex = 0.9, cat.cex = 0.8)

nPOMV24_UP <- length(unique_POMV24$Trend[unique_POMV24$Trend == "UP"]) + length(intersection_24$Trend[intersection_24$Trend == "UP"])
nISAV24_UP <- length(unique_ISAV24$Trend[unique_ISAV24$Trend == "UP"]) + length(intersection_24$Trend[intersection_24$Trend == "UP"])
nINT24_UP <- length(intersection_24$Trend[intersection_24$Trend == "UP"])

UP24 <- draw.pairwise.venn(nPOMV24_UP, nISAV24_UP, nINT24_UP, category = c("POMV 24 hpi", "ISAV 24 hpi"), 
                   lty = rep("blank", 2), fill = c("red", "orange"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = c(0.006, 0.006), scaled = TRUE, 
                   cex = 0.9, cat.cex = 0.8)

nPOMV6_DOWN <- length(unique_POMV6$Trend[unique_POMV6$Trend == "DOWN"]) + length(intersection_6$Trend[intersection_6$Trend == "DOWN"])
nISAV6_DOWN <- length(unique_ISAV6$Trend[unique_ISAV6$Trend == "DOWN"]) + length(intersection_6$Trend[intersection_6$Trend == "DOWN"])
nINT6_DOWN <- length(intersection_6$Trend[intersection_6$Trend == "DOWN"])

DOWN6 <- draw.pairwise.venn(nPOMV6_DOWN, nISAV6_DOWN, nINT6_DOWN, category = c("POMV 6 hpi", "ISAV 6 hpi"), 
                            lty = rep("blank", 2), fill = c("blue", "light blue"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = c(0.006, 0.006), scaled = TRUE, 
                            cex = 0.9, cat.cex = 0.8)

nPOMV24_DOWN <- length(unique_POMV24$Trend[unique_POMV24$Trend == "DOWN"]) + length(intersection_24$Trend[intersection_24$Trend == "DOWN"])
nISAV24_DOWN <- length(unique_ISAV24$Trend[unique_ISAV24$Trend == "DOWN"]) + length(intersection_24$Trend[intersection_24$Trend == "DOWN"])
nINT24_DOWN <- length(intersection_24$Trend[intersection_24$Trend == "DOWN"])

DOWN24 <- draw.pairwise.venn(nPOMV24_DOWN, nISAV24_DOWN, nINT24_DOWN, category = c("POMV 24 hpi", "ISAV 24 hpi"), 
                            lty = rep("blank", 2), fill = c("blue", "light blue"), alpha = rep(0.5, 2), cat.pos = c(0, 180), cat.dist = c(0.006, 0.006), scaled = TRUE, ext.percent = rep(0.01, 3),
                            cex = 0.9, cat.cex = 0.8)

## COMBINE AND SAVE TO FILE
library(gridExtra)
grid.arrange(gTree(children=UP6), gTree(children=DOWN6),
             gTree(children=UP24), gTree(children=DOWN24), ncol=2)


## Get Gene Annotations ##
## get an OrgDb for Atlantic salmon
library(Category)
library(AnnotationHub)
hub <- AnnotationHub()
query(hub, c("salmo salar","orgdb"))
salmodb <- hub[["AH61820"]]
DatPkgFactory(salmodb)
columns(salmodb)

intersection_6_names <- as.character(intersection_6$ENTREZID)
intersection_6_ann <-  AnnotationDbi::select(salmodb, keys = intersection_6_names, columns=c("SYMBOL","GENENAME"), keytype = "ENTREZID")
intersection_6_ann <- merge(intersection_6_ann, intersection_6)
#write.table(intersection_6_ann, 'Results/intersection_6_ann.csv', sep = ',', row.names = FALSE)

intersection_24_names <- as.character(intersection_24$ENTREZID)
intersection_24_ann <-  AnnotationDbi::select(salmodb, keys = intersection_24_names, columns=c("SYMBOL","GENENAME"), keytype = "ENTREZID")
intersection_24_ann <- merge(intersection_24_ann, intersection_24)
#write.table(intersection_24_ann, 'Results/intersection_24_ann.csv', sep = ',', row.names = FALSE)

unique_POMV6_names <- as.character(unique_POMV6$ENTREZID)
unique_POMV6_ann <-  AnnotationDbi::select(salmodb, keys = unique_POMV6_names, columns=c("SYMBOL","GENENAME"), keytype = "ENTREZID")
unique_POMV6_ann <- merge(unique_POMV6_ann, unique_POMV6)
#write.table(unique_POMV6_ann, 'Results/unique_POMV6_ann.csv', sep = ',', row.names = FALSE)

unique_ISAV6_names <- as.character(unique_ISAV6$ENTREZID)
unique_ISAV6_ann <-  AnnotationDbi::select(salmodb, keys = unique_ISAV6_names, columns=c("SYMBOL","GENENAME"), keytype = "ENTREZID")
unique_ISAV6_ann<- merge(unique_ISAV6_ann, unique_ISAV6)
#write.table(unique_ISAV6_ann, 'Results/unique_ISAV6_ann.csv', sep = ',', row.names = FALSE)

unique_POMV24_names <- as.character(unique_POMV24$ENTREZID)
unique_POMV24_ann <-  AnnotationDbi::select(salmodb, keys = unique_POMV24_names, columns=c("SYMBOL","GENENAME"), keytype = "ENTREZID")
unique_POMV24_ann <- merge(unique_POMV24_ann, unique_POMV24)
#write.table(unique_POMV24_ann, 'Results/unique_POMV24_ann.csv', sep = ',', row.names = FALSE)

unique_ISAV24_names <- as.character(unique_ISAV24$ENTREZID)
unique_ISAV24_ann <-  AnnotationDbi::select(salmodb, keys = unique_ISAV24_names, columns=c("SYMBOL","GENENAME"), keytype = "ENTREZID")
unique_ISAV24_ann <- merge(unique_ISAV24_ann, unique_ISAV24)
#write.table(unique_ISAV24_ann, 'Results/unique_ISAV24_ann.csv', sep = ',', row.names = FALSE)


