## Differential Expression Analysis
#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")
library(edgeR)
library(RColorBrewer)
library(pheatmap)
library(statmod)
library(mixOmics)
library(gplots)
library(genefilter)

#setwd("~/Desktop/Scripts")

# import data
countdata <- read.csv('../Data/read_counts_clean.txt', header = TRUE, sep = '\t',  row.names = 'Geneid')

colnames(countdata) <- c('ISAV24.R1', 'ISAV24.R2', 'ISAV24.R3', 'ISAV6.R1', 'ISAV6.R2', 'ISAV6.R3', 
                      'Control.R1', 'Control.R2', 'Control.R3',
                      'POMV24.R1', 'POMV24.R2', 'POMV24.R3', 'POMV6.R1', 'POMV6.R2', 'POMV6.R3')

countdata2 <- countdata[,c(7,8,9,13,14,15,4,5,6,10,11,12,1,2,3)]
head(countdata2)

DataGroups <- c('Control', 'Control', 'Control', 'POMV6', 'POMV6', 'POMV6', 
            'ISAV6', 'ISAV6', 'ISAV6', 'POMV24', 'POMV24', 'POMV24', 
            'ISAV24', 'ISAV24', 'ISAV24')


# Create DGE object of edgeR
dgList <- DGEList(counts=countdata2,group=factor(DataGroups))

y <- dgList

# Filtering and Normalization
# Filter data to retain genes that are represented at least 1 counts per million (cpm) in at least 2 samples
countsPerMillion <- cpm(y)
countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) >= 2)
y <- y[keep,]
y$samples$lib.size <- colSums(y$counts)

# Normalization using TMM method
# TMM normalization is performed to eliminate composition biases between libraries.
y <- calcNormFactors(y, method="TMM")
y$samples

# The performance of the TMM normalization procedure can be examined using mean-difference
# (MD) plots.
# Ideally, the bulk of genes should be centred at a log-fold change of zero. This indicates that
# any composition bias between libraries has been successfully removed. This quality check should
# be repeated by constructing a MD plot for each sample.
plotMD(cpm(y, log=TRUE), column=10)
abline(h=0, col="red", lty=2, lwd=2)

# DATA EXPLORATION
# The data can be explored by generating multi-dimensional scaling (MDS) plots. 
# This visualizes the differences between the expression profiles of different samples in two dimensions.

## MDS plot ## 
#png("plotmds.png")
color <- as.numeric(y$samples$group)
points <- as.numeric(y$samples$group)
plotMDS(y, method="bcv", col = as.numeric(y$samples$group))
#dev.off()

## Hierarchical clustering with heatmaps ##
logcounts <- cpm(y, log=TRUE, prior.count = 1)
head(logcounts)

# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)

# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:1000]
head(select_var)

# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)

head(highly_variable_lcpm)

## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)

# Plot the heatmap
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
colnames(design) <- levels(y$samples$group)
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

# Testing for differential expression
# Fit the linear model
fit <- lmFit(v)
names(fit)

# Contrast Matrix
cont.matrix  <- makeContrasts(ControlvsPOMV6 = Control - POMV6,
                     ControlvsPOMV24 = Control - POMV24,
                     ControlvsISAV6 = Control - ISAV6,
                     ControlvsISAV24 = Control - ISAV24, levels=design)

fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
dim(fit.cont)

cont.matrix2  <- makeContrasts(ControlvsPOMV6 = Control - POMV6,
                              ControlvsISAV6 = Control - ISAV6, levels=design)

fit.cont2 <- contrasts.fit(fit, cont.matrix2)
fit.cont2 <- eBayes(fit.cont2)
dim(fit.cont2)

# We can use the limma decideTests function to generate
# a quick summary of DE genes for the contrasts.

summa.fit <- decideTests(fit.cont)
summary(summa.fit)

summa.fit2 <- decideTests(fit.cont2)
summary(summa.fit2)

vennDiagram(summa.fit2, include=c("up","down"), circle.col = c('green', 'purple'), counts.col = c('green', 'purple'))









# Estimating the dispersion
y <- estimateDisp(y, design, robust=TRUE)
plotBCV(y)

# Fitting the model
fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)
head(fit$coefficients)

# Differential Expression
con <- makeContrasts(Control - POMV6, levels=design)
qlf <- glmQLFTest(fit, contrast=con)
topTags(qlf)

# Here's a closer look at the individual counts-per-million for the top genes. 
top <- topTags(qlf, n=nrow(qlf$table))
head(top$table)
#cpm(y)[top,]

# The total number of genes signicantly up-regulated or down-regulated at 5% FDR is summarized
# as follows:
summary(decideTests(qlf))
plotMD(qlf)

# We use glmTreat to narrow down the list of DE genes and focus on genes that are more biologically
# meaningful.
#tr <- glmTreat(fit, contrast=con, lfc=log2(1.2))
#summary(decideTests(tr))
#plotMD(tr)

# extract and sort differentially expressed genes
sigDownReg <- top$table[top$table$FDR<0.01,]
sigDownReg <- sigDownReg[order(sigDownReg$logFC),]
head(sigDownReg)

sigUpReg <- sigDownReg[order(sigDownReg$logFC, decreasing=TRUE),]
head(sigUpReg)

# Interpreting the DE analysis results
plotSmear(qlf,de.tags = rownames(top$table)[which(top$table$FDR<0.01)])

volcanoData <- cbind(top$table$logFC, -log10(top$table$FDR))
colnames(volcanoData) <- c("logFC", "negLogPval")
head(volcanoData)
plot(volcanoData, pch=19)



anov <- glmQLFTest(fit, contrast=con)
topall <- topTags(anov)









