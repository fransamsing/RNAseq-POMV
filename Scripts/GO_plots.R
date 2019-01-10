
setwd("/Volumes/sam079/RNAseq-POMV")

library(tidyverse)
library(GOplot)
#https://cran.r-project.org/web/packages/GOplot/vignettes/GOplot_vignette.html

## TERMS data set
terms <- read.csv('Results/ego_POMV6.csv')
head(terms)
terms <- terms[,c(1,2,3,7,9)]
colnames(terms) <- c('category', 'ID', 'term', 'adj_pval', 'genes')
terms$genes <- gsub("/", ", ", terms$genes)

## GENES data set
genedata <- read.csv('Results/ControlvsPOMV6_ALL.csv') # CHECK FOR NAs
genes <- genesdata[,c(8,2)]
colnames(genes) <- c('ID', 'logFC')

circ <- circle_dat(terms, genes)

# Generate a simple barplot
GOBar(subset(circ, category == 'BP'))

# Facet the barplot according to the categories of the terms 
GOBar(circ, display = 'multiple')

# Facet the barplot, add a title and change the colour scale for the z-score
GOBar(circ, display = 'multiple', title = 'Z-score coloured barplot', zsc.col = c('yellow', 'black', 'cyan'))

# Generate the bubble plot with a label threshold of 3
GOBubble(circ, labels = 5)

# Add a title, change the colour of the circles, facet the plot according to the categories and change the label threshold
GOBubble(circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3)

# Colour the background according to the category
GOBubble(circ, title = 'Bubble plot with background colour', display = 'multiple', bg.col = T, labels = 3)

# Generate a circular visualization of the results of gene- annotation enrichment analysis
GOCircle(circ)

# Generate a circular visualization for 10 terms
GOCircle(circ, nsub = 10)

# Since we have a lot of significantly enriched processes we selected some specific ones (EC$process)
process <- arrange(terms, adj_pval)
process <- process[1:10,]$term
process <- as.character(process)

# Generate the matrix with selected processes
selected_genes <- arrange(genes, desc(logFC))
selected_genes <- selected_genes[1:500,]

chord <- chord_dat(data = circ, genes = selected_genes, process = process) # Now it is time to generate the binary matrix
head(chord)

chord_df <- data.frame(chord)
chord_df <- rownames_to_column(chord_df)
chord_names <- chord_df$rowname

genenames <- filter(genedata, ENTREZID %in% chord_names)
genenames_vec <- as.character(genenames$Row.names)

rownames(chord) <- genenames_vec

# Create the plot
GOChord(chord, limit = c(0,5), space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 4)

# First, we use the chord object without logFC column to create the heatmap
GOHeat(chord, nlfc = 1, fill.col = c('red', 'yellow', 'green'))

GOCluster(circ, process, clust.by = 'logFC', term.width = 2)


?GOChord



