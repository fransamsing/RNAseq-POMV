## Make heatmap 
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
setwd("/Volumes/sam079/RNAseq-POMV")

## KEGG PATHWAY HEATMAP: RIG-I 
kegg_pathways <- read.csv('Results/kk_POMV24_summary.csv')
rig_I <- as.character(kegg_pathways$geneID[3])
rig_I <- strsplit(rig_I, split = '/')
rig_I <- as.vector(rig_I[[1]])
rig_I

diff_exp_all <- read.csv('Results/Differential_Expression_ALL.csv')
filter_diff_exp_all <- filter(diff_exp_all, ENTREZID %in% rig_I)

filter_diff_exp_all$GENENAME <- as.character(filter_diff_exp_all$GENENAME)
#filter_diff_exp_all[2,11] <- "RIG-I"

# Heatmap
rnames <- filter_diff_exp_all$GENENAME
mat_data <- data.matrix(filter_diff_exp_all[,2:5])
rownames(mat_data) <- rnames

bk1 <- c(seq(-2,0,by=0.1),0.001)
bk2 <- c(0.002,seq(0.5,10,by=0.1))

my_palette <- c(colorRampPalette(colors = c("dodgerblue4", "lightblue", "whitesmoke"))(n = length(bk1)-1),
                c(colorRampPalette(colors = c("floralwhite", "orange", "tomato", "darkred"))(n = length(bk2)-1)))

pheatmap(mat_data, cellwidth = 20, border_color = NA, 
         labels_col = c('POMV 6 hpi', 'POMV 24 hpi', 'ISAV 6 hpi', 'ISAV 24 hpi'), 
         color = my_palette, display_numbers = T, number_color = 'black', number_format = "%.1f")


library(gplots)
library(squash)
my_palette2 <- colorRampPalette(c("navy", "lightblue", "floralwhite", "orange", "tomato", "darkred"))(n = 1000)
heatmap.2(mat_data,
          #cellnote = mat_data,  # same data set for cell labels
          notecex=0.8,
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          col=my_palette2,       # use on color palette defined earlier
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA",
          #Rowv=FALSE,            # turn off column clustering
          key=TRUE, na.color='grey',cexRow=0.9,cexCol=1)

## GO PATHWAY HEATMAP: Defense response to virus (GO:0051607)
GO_pathways <- read.csv('Results/ego_POMV24.csv')

GO_process <- as.character(GO_pathways$geneID[4]) ## 4 = Defense response to virus 
GO_process <- strsplit(GO_process, split = '/')
GO_process <- as.vector(GO_process[[1]])
GO_process  

diff_exp_all <- read.csv('Results/Differential_Expression_ALL.csv')
filter_diff_exp_all <- filter(diff_exp_all, ENTREZID %in% GO_process)
sorted_diff_exp_all <- arrange(filter_diff_exp_all, desc(F))

sorted_diff_exp_all$GENENAME <- as.character(sorted_diff_exp_all$GENENAME)

sorted_diff_exp_all[3,11] <- "interferon-induced GTP-binding protein Mx"

# Heatmap
rnames <- sorted_diff_exp_all$GENENAME
mat_data <- data.matrix(sorted_diff_exp_all[,2:5])
rownames(mat_data) <- rnames

bk1 <- c(seq(-1,-0.5,by=0.1))
bk2 <- c(seq(-0.4,10,by=0.1))

my_palette <- c(colorRampPalette(colors = c("dodgerblue4", "lightblue", "whitesmoke"))(n = length(bk1)-1),
                c(colorRampPalette(colors = c("floralwhite", "orange", "tomato", "darkred"))(n = length(bk2)-1)))

pheatmap(mat_data, cellwidth = 20, border_color = NA, 
         labels_col = c('POMV 6 hpi', 'POMV 24 hpi', 'ISAV 6 hpi', 'ISAV 24 hpi'), 
         color = my_palette, display_numbers = T, number_color = 'black', number_format = "%.1f")


## GO PATHWAY HEATMAP: Immune Response (GO:0006955) 
GO_process <- as.character(GO_pathways$geneID[3]) ## 3 = Immune Response
GO_process <- strsplit(GO_process, split = '/')
GO_process <- as.vector(GO_process[[1]])
GO_process  

filter_diff_exp_all <- filter(diff_exp_all, ENTREZID %in% GO_process)
sorted_diff_exp_all <- arrange(filter_diff_exp_all, desc(F))
sorted_diff_exp_all <- sorted_diff_exp_all[1:30,]
sorted_diff_exp_all$GENENAME <- as.character(sorted_diff_exp_all$GENENAME)

# Heatmap
rnames <- sorted_diff_exp_all$GENENAME
mat_data <- data.matrix(sorted_diff_exp_all[,2:5])
rownames(mat_data) <- rnames

bk1 <- c(seq(-2,0,by=0.1),0.001)
bk2 <- c(0.002,seq(0.5,10,by=0.1))

my_palette <- c(colorRampPalette(colors = c("dodgerblue4", "lightblue", "whitesmoke"))(n = length(bk1)-1),
                c(colorRampPalette(colors = c("floralwhite", "orange", "tomato", "darkred"))(n = length(bk2)-1)))

pheatmap(mat_data, cellwidth = 20, border_color = NA, 
         labels_col = c('POMV 6 hpi', 'POMV 24 hpi', 'ISAV 6 hpi', 'ISAV 24 hpi'), 
         color = my_palette, display_numbers = T, number_color = 'black', number_format = "%.1f")








