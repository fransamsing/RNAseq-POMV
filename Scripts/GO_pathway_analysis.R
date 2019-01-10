# KEGG pathway enrichment analysis
#https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html

#source("https://bioconductor.org/biocLite.R")
#biocLite("clusterProfiler")
biocLite("doseplot")

# Load library and search for the organism
library(clusterProfiler)
library(pathview)
library(tidyverse)
library(enrichplot)

setwd("/Volumes/sam079/RNAseq-POMV")

library(Category)
library(AnnotationHub)
hub <- AnnotationHub()
query(hub, c("salmo salar","orgdb"))
salmodb <- hub[["AH61820"]]
DatPkgFactory(salmodb)
columns(salmodb)

gene_universe <- read.csv('Results/gene_universe.csv')

######################## PATHWAY ANALYSIS FOR POMV 6 HOURS AFTER INFECTION ############################
# Loaded gene lists como from Dif_expression (ALL differentially expressed genes with an adjusted p-value < 0.05)
POMV6 <- read.csv('Results/ControlvsPOMV6_ALL.csv') # CHECK FOR NAs
head(POMV6)
POMV6_SigGeneList <- POMV6[abs(POMV6$logFC) >= 2,]

rownames(POMV6_SigGeneList) <- NULL
head(POMV6_SigGeneList)
POMV6_SigGeneList_logFC <- POMV6_SigGeneList$logFC
names(POMV6_SigGeneList_logFC) <- POMV6_SigGeneList$ENTREZID

ego_POMV6 <- enrichGO(gene = names(POMV6_SigGeneList_logFC),
                universe = as.character(gene_universe$ENTREZID),
                keyType = 'ENTREZID',
                OrgDb         = salmodb,
                ont           = "BP",
                pAdjustMethod = "fdr",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)

summary(ego_POMV6)

barplot(ego_POMV6)
dotplot(ego_POMV6, showCategory=8)
emapplot(ego_POMV6)
cnetplot(ego_POMV6, categorySize="genNUM", foldChange=POMV6_SigGeneList_logFC)

