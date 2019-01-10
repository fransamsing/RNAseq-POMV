# KEGG pathway enrichment analysis
#https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html

#source("https://bioconductor.org/biocLite.R")
#biocLite("clusterProfiler")

# Load library and search for the organism
library(clusterProfiler)
library(pathview)
library(tidyverse)
search_kegg_organism('sasa', by='kegg_code')

######################## PATHWAY ANALYSIS FOR POMV 6 HOURS AFTER INFECTION ############################
# Loaded gene lists como from Dif_expression (ALL differentially expressed genes with an adjusted p-value < 0.05)
POMV6 <- read.csv('Results/ControlvsPOMV6_ALL.csv') # CHECK FOR NAs
head(POMV6)
POMV6_SigGeneList <- POMV6[abs(POMV6$logFC) >= 2,]
POMV6_SigGeneList <- POMV6_SigGeneList[order(abs(POMV6_SigGeneList$logFC), decreasing = TRUE),]
rownames(POMV6_SigGeneList) <- NULL
head(POMV6_SigGeneList)
POMV6_SigGeneList_logFC <- POMV6_SigGeneList$logFC
names(POMV6_SigGeneList_logFC) <- POMV6_SigGeneList$ENTREZID
kk_POMV6 <- enrichKEGG(gene = names(POMV6_SigGeneList_logFC), organism = 'sasa')
head(kk_POMV6)

barplot(kk_POMV6)
dotplot(kk_POMV6)
emapplot(kk_POMV6)
cnetplot(kk_POMV6, categorySize="genNUM", foldChange=POMV6_SigGeneList_logFC)

pathview(gene.data = POMV6_SigGeneList_logFC, 
         pathway.id = "sasa04623", 
         species = "sasa", limit = list(gene = c(-5,5)))

browseKEGG(kk_POMV6, 'sasa04623')

######################## PATHWAY ANALYSIS FOR POMV 24 AFTER INFECTION ############################

# Loaded gene lists como from Dif_expression (ALL differentially expressed genes with an adjusted p-value < 0.05)
POMV24 <- read.csv('Results/ControlvsPOMV24_ALL.csv') # CHECK FOR NAs
head(POMV24)
POMV24_SigGeneList <- POMV24[abs(POMV24$logFC) >= 2,]
POMV24_SigGeneList <- POMV24_SigGeneList[order(abs(POMV24_SigGeneList$logFC), decreasing = TRUE),]
rownames(POMV24_SigGeneList) <- NULL
POMV24_SigGeneList_logFC <- POMV24_SigGeneList$logFC
names(POMV24_SigGeneList_logFC) <- POMV24_SigGeneList$ENTREZID

kk_POMV24 <- enrichKEGG(gene = names(POMV24_SigGeneList_logFC), organism = 'sasa')
head(kk_POMV24)

barplot(kk_POMV24)
dotplot(kk_POMV24)
emapplot(kk_POMV24)
cnetplot(kk_POMV24, categorySize="genNUM", foldChange=POMV24_SigGeneList_logFC)

pathview(gene.data = POMV24_SigGeneList_logFC, 
         pathway.id = "sasa04110", 
         species = "sasa", limit = list(gene = c(-5,5)))

browseKEGG(kk_POMV24, 'sasa04110')

#######################################################################################################
######################## PATHWAY ANALYSIS FOR ISAV 6 HOURS AFTER INFECTION ############################
# Loaded gene lists como from Dif_expression (ALL differentially expressed genes with an adjusted p-value < 0.05)
ISAV6 <- read.csv('Results/RNA_seq/ControlvsISAV6_ALL.csv') # CHECK FOR NAs
head(ISAV6)
ISAV6_SigGeneList <- ISAV6[abs(ISAV6$logFC) >= 2,]
ISAV6_SigGeneList <- ISAV6_SigGeneList[order(abs(ISAV6_SigGeneList$logFC), decreasing = TRUE),]
rownames(ISAV6_SigGeneList) <- NULL
head(ISAV6_SigGeneList)
ISAV6_SigGeneList_logFC <- ISAV6_SigGeneList$logFC
names(ISAV6_SigGeneList_logFC) <- ISAV6_SigGeneList$ENTREZID

kk_ISAV6 <- enrichKEGG(gene = names(ISAV6_SigGeneList_logFC), organism = 'sasa')
head(kk_ISAV6)
barplot(kk_ISAV6)
dotplot(kk_ISAV6)
emapplot(kk_ISAV6)
cnetplot(kk_ISAV6, categorySize="genNUM", foldChange=ISAV6_SigGeneList_logFC)
heatplot(kk_ISAV6) + ggplot2::coord_flip()

pathview(gene.data = ISAV6_SigGeneList_logFC, 
         pathway.id = "sasa04217", 
         species = "sasa", limit = list(gene = c(-5,5)))

######################## PATHWAY ANALYSIS FOR POMV 24 AFTER INFECTION ############################

# Loaded gene lists como from Dif_expression (ALL differentially expressed genes with an adjusted p-value < 0.05)
ISAV24 <- read.csv('Results/ControlvsISAV24_ALL.csv') # CHECK FOR NAs
head(ISAV24)
ISAV24_SigGeneList <- ISAV24[abs(ISAV24$logFC) >= 2,]
ISAV24_SigGeneList <- ISAV24_SigGeneList[order(abs(ISAV24_SigGeneList$logFC), decreasing = TRUE),]
rownames(ISAV24_SigGeneList) <- NULL
head(ISAV24_SigGeneList)
ISAV24_SigGeneList_logFC <- ISAV24_SigGeneList$logFC
names(ISAV24_SigGeneList_logFC) <- ISAV24_SigGeneList$ENTREZID

kk_ISAV24 <- enrichKEGG(gene = names(ISAV24_SigGeneList_logFC), organism = 'sasa')
head(kk_ISAV24)
barplot(kk_ISAV24)
dotplot(kk_ISAV24)
emapplot(kk_ISAV24)
cnetplot(kk_ISAV24, categorySize="genNUM", foldChange=ISAV24_SigGeneList_logFC)

pathview(gene.data = ISAV24_SigGeneList_logFC, 
         pathway.id = "sasa04217", 
         species = "sasa", limit = list(gene = c(-5,5)))


##### NETWORK ANALYSIS USING KEGG PATHWAY ##### 
# Manually download the RIG I pathway from https://www.kegg.jp/kegg-bin/download?entry=sasa04622&format=kgml

RIG_I <- keggGet(dbentries = 'sasa04622')
RIG_I[[1]]$GENE


mapkG <- parseKGML2Graph(RIG_I,expandGenes=TRUE)














