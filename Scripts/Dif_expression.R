## Differential Expression Analysis
#source("https://bioconductor.org/biocLite.R")
#biocLite("Glimma")
library(edgeR)
library(limma)
library(Glimma)
library(RColorBrewer)
library(mixOmics)
library(gplots)
library(genefilter)
library(ggplot2)
library(knitr)
library(pheatmap)
library(reshape2)
library(ggrepel)
library(DEGreport)
library(DESeq2)
library(pheatmap)

# import data
countdata <- read.csv('Data/read_counts_clean.txt', header = TRUE, sep = '\t',  row.names = 'Geneid')
colnames(countdata)

colnames(countdata) <- c('ISAV24.R1', 'ISAV24.R2', 'ISAV24.R3', 'ISAV6.R1', 'ISAV6.R2', 'ISAV6.R3', 
                      'Control.R1', 'Control.R2', 'Control.R3',
                      'POMV24.R1', 'POMV24.R2', 'POMV24.R3', 'POMV6.R1', 'POMV6.R2', 'POMV6.R3')

tail(countdata)
countdata2 <- countdata[,c(7,8,9,13,14,15,4,5,6,10,11,12,1,2,3)]
colnames(countdata2)

DataGroups <- c('Control', 'Control', 'Control', 'POMV6', 'POMV6', 'POMV6', 
            'ISAV6', 'ISAV6', 'ISAV6', 'POMV24', 'POMV24', 'POMV24', 
            'ISAV24', 'ISAV24', 'ISAV24')


# Create DGE object of edgeR
dgList <- DGEList(counts=countdata2,group=factor(DataGroups))
y <- dgList

# Visualizing Library Sizes 
library.size <- as.data.frame(y$samples$lib.size)
library.size$samples <- colnames(y)
colnames(library.size) <- c('library_sizes', 'samples')

ggplot(library.size, aes(samples, library_sizes)) + geom_col() +
  labs(x = "Samples", 
       y= "Library size",
       title = "Barplot of library sizes") + 
  theme(axis.text.x = element_text(face="bold", size=10, angle=90))

# FILTERING DATA
# Filter data to retain genes that are represented at least 0.5 counts per million (cpm) in at least 3 samples
countsPerMillion <- cpm(y)
head(countsPerMillion) # As a general rule, a good threshold can be chosen by identifying the CPM that corresponds to a count of 10, which in this case is about 0.5
countCheck <- countsPerMillion > 0.5
head(countCheck) ## produces a matrix of Trues and False

# Summary of how many TRUEs there are in each row
# There are 11433 genes that have TRUEs in all 12 samples.
table(rowSums(countCheck))

# Let's have a look and see whether our threshold of 0.5 does indeed correspond to a count of about 10-15
# We will look at the first sample
plot(countsPerMillion[,1],countdata2[,1], ylim=c(0,50),xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.5, col = 'blue')
abline(h=10, col = 'blue')
dev.off()

#Apply filter
keep <- which(rowSums(countCheck) >= 3) # 3 samples 
y <- y[keep,]
y$samples$lib.size <- colSums(y$counts)

# NORMALIZATION
# Get log2 counts per million
logcounts <- cpm(y,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

# Normalization using TMM method
# TMM normalization is performed to eliminate composition biases between libraries.
y <- calcNormFactors(y, method="TMM")
y$samples

# The performance of the TMM normalization procedure can be examined using mean-difference
# (MD) plots.
# Ideally, the bulk of genes should be centred at a log-fold change of zero. This indicates that
# any composition bias between libraries has been successfully removed. 
# This quality check should be repeated by constructing a MD plot for each sample.
par(mfrow=c(1,3))
plotMD(cpm(y, log=TRUE), column=10) # the column in indicates what sample I want to make my plot for
abline(h=0, col="red", lty=2, lwd=2)
plotMD(cpm(y, log=TRUE), column=11) 
abline(h=0, col="red", lty=2, lwd=2)
plotMD(cpm(y, log=TRUE), column=10) 
abline(h=0, col="red", lty=2, lwd=2)

dev.off()

# DATA EXPLORATION
# The data can be explored by generating multi-dimensional scaling (MDS) plots. 
# This visualizes the differences between the expression profiles of different samples in two dimensions.

## MDS plot ## 
#png("plotmds.png")
color <- as.numeric(y$samples$group)
points <- as.numeric(y$samples$group)
plotMDS(y, method="bcv", col = as.numeric(y$samples$group))
dev.off()

## Hierarchical clustering with heatmaps ##
logcounts <- cpm(y, log=TRUE, prior.count = 1) # prior.count = average count to be added to each observation to avoid taking log of zero. Used only if log=TRUE.
head(logcounts)

# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes, 100)

# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:100]
head(select_var)

# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)

head(highly_variable_lcpm)

## Get some nicer colours
mypalette <- brewer.pal(11,"PiYG")
morecols <- colorRampPalette(mypalette)

# Plot the heatmap
par(oma=c(2,2,2,2))
colnames(highly_variable_lcpm)
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none",scale="row")

## DIFFERENTIAL EXPRESSION WITH LIMMA-VOOM
# THE DESIGN MATRIX 
# Specify factors
#treatment <-as.factor(c(rep("Control",3),rep("POMV",3),rep("ISAV",3), rep("POMV",3),rep("ISAV",3)))
#time <- as.factor(c(rep("0",3),rep("6",6),rep("24",6)))
#data.frame(Sample=colnames(y),treatment,time)

groups <- y$samples$group
groups

design <- model.matrix(~ 0 + groups)
design
colnames(design)
levels(y$samples$group)
colnames(design) <- levels(y$samples$group)
design
rownames(design)
colnames(y)
rownames(design) <- colnames(y)
design

# Voom transform the data
par(mfrow=c(1,1))
v <- voom(y,design,plot = TRUE)

# What is contained in this object?
names(v)
v$design
v$targets
head(v$E)
v$weights

# The expression values in v$E are already log2 values so we don’t need to log-transform.
par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")

## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")

boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(v$E),col="blue")
dev.off()

# Testing for differential expression
# Fit the linear model
fit <- lmFit(v, method = "robust")
names(fit)

# Contrast Matrix
cont.matrix  <- makeContrasts(ControlvsPOMV6 = POMV6 - Control,
                     ControlvsPOMV24 = POMV24 - Control,
                     ControlvsISAV6 = ISAV6 - Control,
                     ControlvsISAV24 = ISAV24 - Control, levels=design)

fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
dim(fit.cont)

#cont.matrix2  <- makeContrasts(ControlvsPOMV6 = Control - POMV6, ControlvsISAV6 = Control - ISAV6, levels=design)
#fit.cont2 <- contrasts.fit(fit, cont.matrix2)
#fit.cont2 <- eBayes(fit.cont2)
#dim(fit.cont2)

# We can use the limma decideTests function to generate
# a quick summary of DE genes for the contrasts.

summa.fit <- decideTests(fit.cont)
summary(summa.fit)

expres <- summary(summa.fit)
expres <- t(expres)
expres <- as.data.frame(expres)
expres
myColors <- c("#A1D76A", "#FFFF99", "#E9A3C9")
ggplot(expres, aes(x = Var1, y = Freq)) + geom_col(aes(fill = Var2)) +  
  labs(x = "Comparison between samples", 
  y= "Number of genes",
  color = "Change in gene expression",
  title = "Up or down regulation of my genes") +
  scale_fill_manual(values=myColors)

#summa.fit2 <- decideTests(fit.cont2)
#summary(summa.fit2)
#vennDiagram(summa.fit2, include=c("up","down"), circle.col = c('green', 'purple'), counts.col = c('green', 'purple'))

# The toptable
## To get the full table
#ControlvsPOMV6_fulltable <- topTable(fit.cont,coef="ControlvsPOMV6",sort.by="p",n="Inf")

sig.genes <- function(fit.contrast, mycoef, upregulated = c('UP', 'DOWN', 'ALL')) {
  x <- topTable(fit = fit.contrast, coef = mycoef, sort.by="p",n="Inf")
  x_sig <- x[x$adj.P.Val < 0.05,] 
  if (upregulated == 'UP') {
    x_UP <- x_sig[x_sig$logFC >= 1.5,]
    #x_UP <- rownames_to_column(x_UP, var = "GeneID")
  } else if (upregulated == 'DOWN') {
    x_DOWN <- x_sig[x_sig$logFC <= -1.5,]
    #x_DOWN <- rownames_to_column(x_DOWN, var = "GeneID")
  }
  else {
   x_sig 
  }
  }

# POMV
ControlvsPOMV6_UP <- sig.genes(fit.cont, mycoef = "ControlvsPOMV6", 'UP')
ControlvsPOMV6_DOWN <- sig.genes(fit.cont, mycoef = "ControlvsPOMV6", 'DOWN')
ControlvsPOMV6_ALL <- sig.genes(fit.cont, mycoef = "ControlvsPOMV6", upregulated = FALSE)

ControlvsPOMV24_UP <- sig.genes(fit.cont, mycoef = "ControlvsPOMV24", 'UP')
head(ControlvsPOMV24_UP[order(-ControlvsPOMV24_UP$logFC),])
ControlvsPOMV24_DOWN <- sig.genes(fit.cont, mycoef = "ControlvsPOMV24", 'DOWN')
ControlvsPOMV24_ALL <- sig.genes(fit.cont, mycoef = "ControlvsPOMV24", upregulated = FALSE)

## ISAV
ControlvsISAV6_UP <- sig.genes(fit.cont, mycoef = "ControlvsISAV6", 'UP')
ControlvsISAV6_DOWN <- sig.genes(fit.cont, mycoef = "ControlvsISAV6", 'DOWN')
ControlvsISAV6_ALL <- sig.genes(fit.cont, mycoef = "ControlvsISAV6", upregulated = FALSE)

ControlvsISAV24_UP <- sig.genes(fit.cont, mycoef = "ControlvsISAV24", 'UP')
ControlvsISAV24_DOWN <- sig.genes(fit.cont, mycoef = "ControlvsISAV24", 'DOWN')
ControlvsISAV24_ALL <- sig.genes(fit.cont, mycoef = "ControlvsISAV24", upregulated = FALSE)

# PLOTS FOR DE genes

# THE VOLCANO PLOT
res_table <- topTable(fit = fit.cont, coef = 1, sort.by="p",n="Inf")
## Obtain logical vector regarding whether padj values are less than 0.05
threshold <- res_table$adj.P.Val < 0.05 
## Determine the number of TRUE values
length(which(threshold))
## Add logical vector as a column (threshold) to the res_tableOE
res_table$threshold <- threshold

## Volcano plot
ggplot(res_table) +
  geom_point(aes(x=logFC, y=-log10(adj.P.Val), colour=threshold)) +
  ggtitle("Control vs POMV6") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_y_continuous(limits = c(0,25)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 

# Add some labels
## Sort by ordered padj
#res_table_ordered <- res_table[order(res_table$adj.P.Val), ] 
## Create a column to indicate which genes to label
#res_table_ordered$genelabels <- ""
#res_table_ordered$genelabels[1:15] <- rownames(res_table_ordered)[1:15]
#ggplot(res_table_ordered) +
  #geom_point(aes(x=logFC, y=-log10(adj.P.Val), colour=threshold)) +
  #geom_text_repel(aes(x = logFC, y = -log10(adj.P.Val), label = ifelse(genelabels == T, rownames(res_table_ordered),""))) +
  #ggtitle("Control vs POMV") +
  #xlab("log2 fold change") + 
  #ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,30)) +
  #scale_x_continuous(limits = c(-10, 10)) +
  #theme(legend.position = "none",
        #plot.title = element_text(size = rel(1.5), hjust = 0.5),
        #axis.title = element_text(size = rel(1.25)))

## MA Plots
head(fit.cont$coefficients)
limma::plotMA(fit.cont,coef=1,status=summa.fit[,"ControlvsPOMV6"], values = c(-1, 1), hl.col=c("blue","red"), hl.pch = 16, bg.pch = 16, bg.cex = 0.3)

## MA Plots in Glimma
#summary(summa.fit)
#glMDPlot(fit.cont, status=summa.fit, counts=v, groups=groups, side.main="Symbols")

## Gene Ontology Analysis ##
## get an OrgDb for Atlantic salmon
library(Category)
library(AnnotationHub)
hub <- AnnotationHub()
query(hub, c("salmo salar","orgdb"))
salmodb <- hub[["AH61820"]]
DatPkgFactory(salmodb)
columns(salmodb)


#https://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf

ann <-  AnnotationDbi::select(salmodb, keys=rownames(fit.cont), columns=c("ENTREZID","GENENAME"), keytype = "SYMBOL")
head(ann, 10)

# control vs 6 hours - POMV
ControlvsPOMV6_UP_ann <- merge(ControlvsPOMV6_UP, ann, by.x=0, by.y="SYMBOL")
head(ControlvsPOMV6_UP_ann)

ControlvsPOMV6_DOWN_ann <- merge(ControlvsPOMV6_DOWN, ann, by.x=0, by.y="SYMBOL")
head(ControlvsPOMV6_DOWN_ann)

ControlvsPOMV6_ALL_ann <- merge(ControlvsPOMV6_ALL, ann, by.x=0, by.y="SYMBOL")
head(ControlvsPOMV6_ALL_ann)

## control vs 24 hours - POMV
ControlvsPOMV24_UP_ann <- merge(ControlvsPOMV24_UP, ann, by.x=0, by.y="SYMBOL")
#head(ControlvsPOMV24_UP_ann)
head(ControlvsPOMV24_UP_ann[order(-ControlvsPOMV24_UP_ann$logFC),])

ControlvsPOMV24_DOWN_ann <- merge(ControlvsPOMV24_DOWN, ann, by.x=0, by.y="SYMBOL")
head(ControlvsPOMV24_DOWN_ann)

ControlvsPOMV24_ALL_ann <- merge(ControlvsPOMV24_ALL, ann, by.x=0, by.y="SYMBOL")
head(ControlvsPOMV24_ALL_ann)

# control vs 6 hours - ISAV
ControlvsISAV6_UP_ann <- merge(ControlvsISAV6_UP, ann, by.x=0, by.y="SYMBOL")
head(ControlvsISAV6_UP_ann)

ControlvsISAV6_DOWN_ann <- merge(ControlvsISAV6_DOWN, ann, by.x=0, by.y="SYMBOL")
head(ControlvsISAV6_DOWN_ann)

ControlvsISAV6_ALL_ann <- merge(ControlvsISAV6_ALL, ann, by.x=0, by.y="SYMBOL")
head(ControlvsISAV6_ALL_ann)

# control vs 24 hours - ISAV
ControlvsISAV24_UP_ann <- merge(ControlvsISAV24_UP, ann, by.x=0, by.y="SYMBOL")
head(ControlvsISAV24_UP_ann)

ControlvsISAV24_DOWN_ann <- merge(ControlvsISAV24_DOWN, ann, by.x=0, by.y="SYMBOL")
head(ControlvsISAV24_DOWN_ann)

ControlvsISAV24_ALL_ann <- merge(ControlvsISAV24_ALL, ann, by.x=0, by.y="SYMBOL")
head(ControlvsISAV24_ALL_ann)

# Write tables to results
write.table(ControlvsPOMV6_UP_ann, 'Results/RNA_seq/ControlvsPOMV6_UP.csv', row.names = FALSE, sep = ",")
write.table(ControlvsPOMV6_DOWN_ann, 'Results//RNA_seq/ControlvsPOMV6_DOWN.csv', row.names = FALSE, sep = ",")
write.table(ControlvsPOMV6_ALL_ann, 'Results//RNA_seq/ControlvsPOMV6_ALL.csv', row.names = FALSE, sep = ",")

write.table(ControlvsPOMV24_UP_ann, 'Results/RNA_seq/ControlvsPOMV24_UP.csv', row.names = FALSE, sep = ",")
write.table(ControlvsPOMV24_DOWN_ann, 'Results/RNA_seq/ControlvsPOMV24_DOWN.csv', row.names = FALSE, sep = ",")
write.table(ControlvsPOMV24_ALL_ann, 'Results/RNA_seq/ControlvsPOMV24_ALL.csv', row.names = FALSE, sep = ",")

write.table(ControlvsISAV6_UP_ann, 'Results/RNA_seq/ControlvsISAV6_UP.csv', row.names = FALSE, sep = ",")
write.table(ControlvsISAV6_DOWN_ann, 'Results/RNA_seq/ControlvsISAV6_DOWN.csv', row.names = FALSE, sep = ",")
write.table(ControlvsISAV6_ALL_ann, 'Results/RNA_seq/ControlvsISAV6_ALL.csv', row.names = FALSE, sep = ",")

write.table(ControlvsISAV24_UP_ann, 'Results/RNA_seq/ControlvsISAV24_UP.csv', row.names = FALSE, sep = ",")
write.table(ControlvsISAV24_DOWN_ann, 'Results/RNA_seq/ControlvsISAV24_DOWN.csv', row.names = FALSE, sep = ",")
write.table(ControlvsISAV24_ALL_ann, 'Results/RNA_seq/ControlvsISAV24_ALL.csv', row.names = FALSE, sep = ",")
