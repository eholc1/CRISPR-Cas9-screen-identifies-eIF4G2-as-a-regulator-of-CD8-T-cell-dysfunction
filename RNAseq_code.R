#Bulk RNA-seq analysis

library(DESeq2)
library(dplyr)
library(magrittr)
library(RColorBrewer)
library(ggrepel)
library(biomaRt)
library(pheatmap)
library(ComplexHeatmap)

#### Data import ####
# Note if using R < 4.0.0, set stringsAsFactors = FALSE in read.delim
# DESeq2 needs un-normalized counts, so don't use TPM or FPKM
data_orig <- read.delim("gene_expected_count.annot.txt", row.names = NULL)
data <- data_orig

#Reorder columns 
data <- data[, c('gene_id', 'entrezgene_id', 'external_gene_name', 'description',
                 'X8430.AH.1', 'X8430.AH.2', 'X8430.AH.3', 
                 'X8430.AH.4', 'X8430.AH.5', 'X8430.AH.6', 
                 'X8430.AH.7', 'X8430.AH.8' ,'X8430.AH.9', 
                 'X8430.AH.10', 'X8430.AH.11', 'X8430.AH.12')] 

#Rename sample columns
data <- rename(data,
               Acute.Cont.1 = X8430.AH.1,
               Acute.Cont.2 = X8430.AH.2,
               Acute.Cont.3 = X8430.AH.3,
               Acute.Eif.1 = X8430.AH.4,
               Acute.Eif.2 = X8430.AH.5,
               Acute.Eif.3 = X8430.AH.6,
               Chronic.Cont.1 = X8430.AH.7,
               Chronic.Cont.2 = X8430.AH.8,
               Chronic.Cont.3 = X8430.AH.9,
               Chronic.Eif.1 = X8430.AH.10,
               Chronic.Eif.2 = X8430.AH.11,
               Chronic.Eif.3 = X8430.AH.12)

# Deal with genes that don't have annotated gene symbols (external_gene_name)
# Use ENSEMBL ID if gene symbol not available
data$external_gene_name <- ifelse(
  data$external_gene_name == ".",
  data$gene_id,
  data$external_gene_name)

# Deal with duplicated gene symbols
# Combine gene symbol with ENSEMBL ID if non-unique
data$external_gene_name <- ifelse(
  duplicated(data$external_gene_name),
  paste(data$external_gene_name, data$gene_id, sep="_"),
  data$external_gene_name)

# Then we can use the gene symbol column as the row names,
# and subset the count data for further analysis
rownames(data) <- data$external_gene_name
count.data <- data[,5:ncol(data)] # All columns after 4 are count data



#### DESeq2 preparation ####
#Trim count matrix file to include just sample columns
data <- data[ ,-c(1,2,3,4)]

#Make Sample Information Table
meta <- data.frame('samples'=colnames(data)) #Deletes informational columns and leaves sample names

#Assign sample names
#Put control group first (can add levels)
meta$condition <- factor(c(rep(c("Acute_Control"), 3), rep(c("Acute_sgEif4g2"), 3), 
                           rep(c("Chronic_Control"), 3), rep(c("Chronic_sgEif4g2"), 3))) 

row.names(meta) <- meta$samples # set the row names to be the sample names
meta$samples <- NULL # clear the sample name column since isnt needed

#Make sure all column names in count matrix match the row names in the sample label file
all(colnames(data) == rownames(meta))
write.csv(data, file="data.csv")



#### DESeq2 input ####
## Create DESeq object, line by line
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = meta,
                              design = ~ condition)
str(dds)

#Filter out any genes that are lacking counts across any of the samples, i.e.: all columns are zero.
keep <- rowSums(counts(dds)) >= 1
dds <- dds[keep,]
#Raw count table
write.csv(counts(dds, normalized = FALSE), file="DESeq2_raw_counts.csv")

#Normalization
rld <- rlog(dds, blind = TRUE)
#Normalized rlog count table
write.csv(assay(rld), file="DESeq2_rlogNormalized_counts.csv")


#DESeq2 Model Fitting
dds <- DESeq(dds)

resultsNames(dds)
# "Intercept"                                   "condition_Acute_sgEif4g2_vs_Acute_Control"  
# "condition_Chronic_Control_vs_Acute_Control"  "condition_Chronic_sgEif4g2_vs_Acute_Control"


#### QC: COUNT BOXPLOTS ####
## Setup for raw counts
pdata = data.frame(colData(dds))
mat = as.matrix(assay(dds))
title = 'Raw counts'
y_label = 'log2(counts)'

# Create annotationn table for raw plots
annot_df = data.frame(
  sample = row.names(pdata),
  Group = factor(pdata[, "condition"]),
  row.names = row.names(pdata),
  stringsAsFactors = F)

# Join counts and annotation table
tidy_mat = tidyr::gather(as_tibble(mat), key = 'sample', value = 'counts') %>%
  left_join(annot_df, by = 'sample')

#Plot
box_plot = ggplot(tidy_mat, aes(x = sample, y = log2(counts), fill = Group)) +
  geom_boxplot(notch = TRUE) +
  labs(
    title = title,
    x = '',
    y = y_label) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90))
box_plot

#To understand how the rlog normalization impacted the distributions of counts for each sample, 
#we can plot boxplots for the normalized data and compare that to our plot of the raw data.
## rlog counts
pdata = data.frame(colData(rld))
mat = as.matrix(assay(rld))
title = 'Rlog normalized counts'
y_label = 'rlog(counts)'

annot_df = data.frame(
  sample = row.names(pdata),
  Group = factor(pdata[, "condition"]),
  row.names = row.names(pdata),
  stringsAsFactors = F
)

tidy_mat = tidyr::gather(as_tibble(mat), key = 'sample', value = 'counts') %>%
  left_join(annot_df, by = 'sample')

group_colors <- list(Group=c(Acute_Control="snow3", Acute_sgEif4g2="lightskyblue", Chronic_Control='black', Chronic_sgEif4g2="blue"))

box_plot = ggplot(tidy_mat, aes(x = sample, y = counts, fill = Group)) +
  scale_fill_manual(values=c(Acute_Control="snow3", Acute_sgEif4g2="lightskyblue", Chronic_Control='black', Chronic_sgEif4g2="blue")) +
  geom_boxplot(notch = TRUE) +
  labs(
    title = title,
    x = '',
    y = y_label) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90))
box_plot



#### DEG  ####
#check what comparisons were automatically generated during fitting
resultsNames(dds)

Comparison <- "condition_Acute_sgEif4g2_vs_Acute_Control"
Comparison <- "condition_Chronic_Control_vs_Acute_Control"
Comparison <- "condition_Chronic_sgEif4g2_vs_Acute_Control"

#Choose comparison groups
res <- results(dds, contrast=c("condition", "Acute_sgEif4g2", "Acute_Control"))
res <- results(dds, contrast=c("condition", "Chronic_Control", "Acute_Control")) 
res <- results(dds, contrast=c("condition", "Chronic_sgEif4g2", "Chronic_Control"))

head(res)


## ASSIGN THRESHOLDS
fc <- 1.5
pval <- 0.05

df<- res[order(res$padj),] #select our data of interest
df <- as.data.frame(df) #convert our object type
df <- cbind("id" = row.names(df), df) #set rownames to valid column
str(df) 

#Subset data to label the datapoints (genes) that pass our thresholds.
df$dot <- rep(3, nrow(df))
df$dot[which(df$padj <= pval & df$log2FoldChange < 0 & abs(df$log2FoldChange) >= log2(fc))] = 2
df$dot[which(df$padj <= pval & df$log2FoldChange > 0 & abs(df$log2FoldChange) >= log2(fc))] = 1
# 1 = UP, 2= DOWN

#Count the number of significant up and down genes
df$dot <- factor(df$dot,levels = c(1,2,3), labels = c(paste0('Up: ', sum(df$dot == 1)),paste0('Down: ', sum(df$dot == 2)),'NS'))

#Write csv to view all DE genes
write.csv(df, file="Chronic.Cont_vs_Acute.Cont2.csv")
write.csv(df, file="Acute.Eif_vs_Acute.Cont.csv")
write.csv(df, file="Chronic.Eif_vs_Chronic.Cont2.csv")



#### HEATMAP (Figure 4D & S4B) ####
color_palette <- colorRampPalette(brewer.pal(9, "RdBu"))(100)
reversed_palette <- rev(color_palette)
group_colors <- list(Group=c(Acute_Control="snow3", Acute_sgEif4g2="lightskyblue", Chronic_Control='black', Chronic_sgEif4g2="blue"))

#Figure S4D: select the top 2000 expressed genes across all samples
select <- order(rowMeans(assay(rld)), decreasing=TRUE)[1:2000]

#Figure 4D: select specific genes 
select <- c('Stmn1', 'Cdk1', 'Cdk2', 'Bub1b', 'Mcm3', 'Mcm4', 'Cdca2', 'Ccnb1', 'E2f2')

select <- c('Slc7a11', 'Aifm2', 'Alox5', 'Ncoa4', 'Cd36', 'Gclc', 'Acsl3', 'Mboat7', 'Prdx4', 'Gpx7', 'Apaf1', 'Pycard', 'Fas')

select <- c('Hspbap1', 'Hspa4l', 'Hspb6', 'Hspb11')

select <- c('Il7r', 'Foxo1', 'Ccr7', 'Id3', 'Il15ra', 'Sell')

select <- c( 'Cd44', 'Il2ra', 'Gzme', 'Gzmg', 'Gzmd', 'Gzmc', 'Ccl3', 'Ccl4', 'Cd69', 'Ifng', 'Tnf', 'Gzmk', 'Prf1', 'Il10', 'Myc', 'Arid1a','Prdm1', 'Klrg1', 'Irf1', 'Irf4', 'Irf6', 'Irf7', 'Irf8', 'Nfkb1', 'Plcg1',
             'Il16', 'Il24', 'Stat1', 'Stat2', 'Stat4', 'Il18r1', 'Il12rb2', 'Icos', 'Tbx21', 'Batf', 'Batf3', 'Cxcr6')

select <- c('Gzmc', 'Gzmd', 'Gzme', 'Gzmg', 'Gzmk', 'Prf1', 'Il16', 'Tnf', 'Ccl3', 'Ccl4', 'Il10', 'Il24', 'Ifng') #Cytokines

select <- c('Cd44', 'Il2ra', 'Il18r1', 'Il12rb2', 'Icos',
            'Plcg1', 'Myc', 'Arid1a', 'Irf4', 
            'Irf8', 'Klrg1', 'Prdm1', 'Batf') #Effector

select <- c('Pdcd1', 'Havcr2', 'Tigit', 'Tnfrsf9', 'Entpd1', 'Cd38', 'Cd244a', 'Cd69',
            'Cd101', 'Cd160', 'Nfatc1', 'Tox', 'Bhlhe40',
            'Ctla4', 'Lag3')

select <- c('Cx3cr1', 'Zeb2',  'Klra1', 'Klra3', 'Klra5', 'Klra6', 'Klra17', 'Klrc1', 'Klrc2', 'Klri2', 'Klrc3', 'Klrd1', 'Klre1', 'Klrh1', 'Klrk1',
            'Klrb1', 'Klrb1a', 'Klrb1b', 'Klrb1c')


#The pheatmap function does quite a lot in a single step, 
#including scaling the data by row and clustering both the samples (columns) and genes (rows).
df <- data.frame(Group = colData(rld)[,c('condition')], row.names = rownames(colData(dds)))

p <- pheatmap(assay(rld)[select,], scale="row",  cluster_rows=TRUE, show_rownames=FALSE, 
              cluster_cols=FALSE, annotation_col=df, color = reversed_palette,
              labels_col=c('Acute.Cont.1', 'Acute.Cont.2', 'Acute.Cont.3', 'Acute.Eif.1', 'Acute.Eif.2', 'Acute.Eif.3',
                           'Chronic.Cont.1', 'Chronic.Cont.2', 'Chronic.Cont.3', 'Chronic.Eif.1', 'Chronic.Eif.2', 'Chronic.Eif.3'),
              annotation_colors=group_colors, fontsize = 9, fontsize_row = 9, main = 'Transcriptomics Top 2000')

p <- pheatmap(assay(rld)[select,], scale="row",  
              cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=TRUE, annotation_col=df, show_row_dend = FALSE, 
              labels_col=c('Acute_Control.1', 'Acute_Control.2', 'Acute_Control.3', 'Acute_sgEif4g2.1', 'Acute_sgEif4g2.2', 'Acute_sgEif4g2.3',
                           'Chronic_Control.1', 'Chronic_Control.2', 'Chronic_Control.3', 'Chronic_sgEif4g2.1', 'Chronic_sgEif4g2.2', 'Chronic_sgEif4g2.3'),
              annotation_colors=group_colors, color = reversed_palette, fontsize = 10, fontsize_row = 10, cellwidth=10, cellheight = 10)
dev.off()
p

#Add vertical lines to heatmap
co=column_order(p)
nc = ncol(p)

decorate_heatmap_body(heatmap = "matrix_1", code = {
  # Define the number of groups
  n_groups = nc %/% 3  # Integer division to get full groups of 3
  for (group in 0:(n_groups - 1)) {
    i_start = group * 3 + 1
    i_end = min((group + 1) * 3, nc)  # Ensure we don't go beyond the number of columns
    # Draw a rectangle around the current group of columns
    grid.rect(x = (i_start + i_end - 1) / 2 / nc,  # Center position
              width = (i_end - i_start + 1) / nc,  # Width covering the group
              height = unit(1, "npc"),  # Full height of the heatmap
              gp = gpar(col = "black", fill = NA, lwd = 0.75))  # Border properties
  }
})

decorate_annotation('Group', code = {
  # Define the number of groups
  n_groups = nc %/% 3  # Integer division to get full groups of 3
  for (group in 0:(n_groups - 1)) {
    i_start = group * 3 + 1
    i_end = min((group + 1) * 3, nc)  # Ensure we don't go beyond the number of columns
    # Draw a rectangle around the current group of columns
    grid.rect(x = (i_start + i_end - 1) / 2 / nc,  # Center position
              width = (i_end - i_start + 1) / nc,  # Width covering the group
              height = unit(1, "npc"),  # Full height of the heatmap
              gp = gpar(col = "black", fill = NA, lwd = 0.75))  # Border properties
  }
})



#### PCA Plot (Figure S4A) ####
#PCA plot for Rlog-Normalized counts for all samples
Group <- factor(meta$condition)
SampleName <- factor(row.names(meta))

#Generate the PCA projections
p.all <- plotPCA(rld, intgroup = c('condition'), ntop = 500)

gp <- ggplot(p.all$data, aes(x = PC1, y = PC2), color = Group, shape = Group) +
  xlab(p.all$labels[2]) + ylab(p.all$labels[1]) +
  geom_point(aes(fill=Group, shape = Group), shape = 21, size=5, color = 'black') + #Shape refers to what kind of pointshape/outline (21 is outlined circle)
  scale_fill_manual(values=c(Acute_Control="snow3", Acute_sgEif4g2="lightskyblue", 
                             Chronic_Control='black', Chronic_sgEif4g2="blue")) +
  ggtitle(label = as.character('Transcriptomics'))  + theme_bw(base_size = 10) +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = c(0.79,0.81), legend.key.size = unit(0.4, "cm"),
        legend.background = element_rect(fill="white", size=0.5, linetype="solid", colour ="black"))

plot(gp)

ggsave('PCA_plot.jpg', plot = gp, width = 4, height = 4)




#### PATHWAY ANALYSIS (Figure 4C & S4F) ####
library(tidyverse)
library(org.Mm.eg.db)
library(ReactomePA)
library(msigdbr)
library(ComplexHeatmap)
library(clusterProfiler)
library(enrichplot)

#Choose comparison
res <- results(dds, contrast=c("condition", "Chronic_Control", "Acute_Control")) 
res <- results(dds, contrast=c("condition", "Acute_sgEif4g2", "Acute_Control"))
res <- results(dds, contrast=c("condition", "Chronic_sgEif4g2", "Chronic_Control"))

res <- res[order(-res$log2FoldChange), ] #Order in decreasing LFC (required for Clusterprofiler)
gene_list <- res$log2FoldChange #Make list with just gene name and LFC
names(gene_list) <- row.names(res)
gene_list<-na.omit(gene_list)


## MSIGDB: Hallmark
gsea_h_df <- data.frame()
m_df <- msigdbr(species = 'Mus musculus', category = 'H') %>% dplyr::select(gs_name, gene_symbol)
gsea_h <- GSEA(gene_list, TERM2GENE = m_df, pAdjustMethod = 'BH', pvalueCutoff = 0.05)
gsea_h_df <- gsea_h@result

#Figure 4C:
gseaplot2(gsea_h, geneSetID = c(1, 2, 3), pvalue_table = FALSE, subplots = 1:2, base_size = 15,
          color = c("palevioletred1", "mediumseagreen", "dodgerblue2")) 

#Figure S4F:
gseaplot2(gsea_h, geneSetID = c(1, 2, 5), pvalue_table = FALSE, subplots = 1:2, base_size = 15,
          color = c("palevioletred1", "mediumseagreen", "dodgerblue2")) 



#### fGSEA (Figure S4E) ####
library(fgsea)
library(data.table)
library(ggplot2)
library(tibble)

res <- results(dds, contrast=c("condition", "Chronic_Control", "Acute_Control")) 

#Make gene set file that contains one column with gene name (id) and one column with a ranked metric (logFC)
res <- na.omit(res)
#res <- res[res$baseMean>50, ]
df <- res[order(-res$log2FoldChange), ] #Order in decreasing LFC
logFC <- df$log2FoldChange #Make list with just gene name and LFC
names(logFC) <- row.names(df)
ranks <- as.data.frame(logFC) #convert our object type
ranks <- tibble::rownames_to_column(ranks, "id") #Change row names to column

#Convert logFC values to list format, then to numeric values
ranks_list <- list(ranks$logFC)
ranks_numeric_list <- as.numeric(unlist(ranks_list))

#Name the rank values with the gene names
names(ranks_numeric_list) <- ranks$id


#Read in Wu/Wherry gene set:
gene_set <- read.csv(file='Wu_Wherry_C7vsA7.csv', header=TRUE, sep=',')
gene_set <- gene_set[gene_set$adj.P.Val < 0.05 & abs(gene_set$logFC) >= log2(1.5), ]
gene_set$logFC <- gene_set$logFC * (-1)
gene_set <- gene_set[,c("gene","logFC")]  
gene_set_pos <- gene_set[gene_set$logFC>0, ]
gene_set_pos <- t(gene_set_pos)
gene_set_pos <- as.character(gene_set_pos)
gene_sets <- list(gene_set_pos)
write_csv(gene_set_pos, file = 'Wu_Wherry_pos.csv')


#Run fGSEA
fgseaRes <- fgsea(pathways = gene_sets, 
                  stats = ranks_numeric_list,
                  minSize=15,
                  maxSize=500)

#Plot GSEA for a select gene set (in this case I only have one)
plotEnrichment(gene_sets[[1]],
               ranks_numeric_list) + labs(title="Tex_signature")

dev.off()



