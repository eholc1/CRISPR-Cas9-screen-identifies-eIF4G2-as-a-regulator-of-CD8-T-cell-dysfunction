library(ggplot2)
library(ggrepel)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(ComplexHeatmap)
library(DEqMS)
library(matrixStats)
library(limma)
library(tibble)

#### DATA PROCESSING ####
##Read protein table as input and filter it
df.prot = read.csv("Proteomics_proteins.csv", header=T)

# Extract columns of normalized counts 
data = df.prot[,c(1:14, 17)]
data <- data %>% replace(is.na(.), 0)
data <- data[data$Gene.Symbol != "", ] #Remove unidentified proteins

#Check for duplicate genes
has_duplicates <- any(duplicated(data$Gene.Symbol))
print(has_duplicates)
duplicate_genes <- data$Gene.Symbol[duplicated(data$Gene.Symbol)]
print(duplicate_genes)

# Deal with duplicated gene symbols by combining gene symbol with ENSEMBL ID if non-unique
data$Gene.Symbol <- ifelse(
  duplicated(data$Gene.Symbol),
  paste(data$Gene.Symbol, data$Gene.ID, sep="_"),
  data$Gene.Symbol)

data_psm <- data
## STOP HERE FOR DEqMS PSM TABLE ##


# Then we can use the gene symbol column as the row names
rownames(data) <- data$Gene.Symbol
#Subset the count data for further analysis
data <- data[,3:14] 
#Remove rows with all 0's
data <- data[rowSums(data != 0, na.rm = TRUE) > 0, ]


#### VENN DIAGRAM: RNA vs Proteins (Figure 4A) ####
#VennDiagram
library(VennDiagram)
library(scales) # Required for alpha()
library(grid)

venn_prot <- read.csv(file = 'protein_counts.csv', header=TRUE, sep=',')
venn_RNA <- read.csv(file = 'DESeq2_raw_counts.csv', header=TRUE, sep=',')

ven_list <- list(RNA = venn_RNA$ID, Protein = venn_prot$ID)

p <- venn.diagram(ven_list, 
                  category.names = c('RNA', 'Protein'),
                  filename = 'RNA & Protein Overlap.png', 
                  output = TRUE, imagetype="png", height = 1000, width = 1000, compression = "lzw",
                  lwd = 0.6, col = c('grey40', 'grey40'), fill = c('indianred1', "springgreen2"), alpha=0.5, #Circles
                  cex = 0.5, fontfamily = "sans", #Numbers
                  cat.cex = 0.4, cat.fontface = "bold", cat.fontfamily = "sans",  #Name appearance
                  cat.default.pos = "outer", cat.pos = c(-20, 30), cat.dist = c(0.1, 0.2), margin = 0.1, #Name position
                  ext.pos = 90, ext.length = 0.5, ext.dist = -0.05, ext.line.lwd = 1, ) # Label line
p

plot(1, 1, type = "n", xlab = "", ylab = "", xlim = c(0, 2), ylim = c(0, 2), axes = FALSE)
legend(
  "center", pch=21, pt.bg = c("indianred1", "springgreen2"),pt.cex=1.5,
  legend = c("25438 Identified Transcripts", "7444 Identified Proteins"), 
)

dev.off()


#### DEqMS input ####
dev.off()
par(mfrow = c(4, 3))  # Layout for 2 rows and 5 columns (adjust as needed)
for (i in 1:ncol(data.log)) {
  hist(data.log[, i],
       main = colnames(data.log)[i],  # Title for each sample
       xlab = "Log Intensity",        # X-axis label
       ylab = "Frequency",            # Y-axis label
       col = "lightblue",             # Color
       border = "black")
}

data <- data %>% replace(is.na(.), 0) #Replace NA's with 0
data.log = log2(data+1) #log transform
data.log = na.omit(data.log) #remove rows with NAs

#Use boxplot to check if the samples have medians centered. if not, do median centering.
boxplot(data.log,las=2,main="TMT data")
data.log = equalMedianNormalization(data.log)
boxplot(data.log,las=2,main="TMT data median normalized")

#Make design table
cond = as.factor(c("AC","AC","AC","AE","AE","AE", "CC","CC","CC","CE","CE","CE"))
design = model.matrix(~0+cond) # 0 means no intercept for the linear model
colnames(design) = gsub("cond","",colnames(design))

#Make contrasts
x <- c("AE-AC", "CC-AC", "CE-CC")
cont =  makeContrasts(contrasts=x,levels=design)
fit1 = lmFit(data.log,design = design)
fit2 = contrasts.fit(fit1,contrasts = cont)
fit3 <- eBayes(fit2)

#Assign a extra variable `count` to fit3 object, telling how many PSMs are quantifed for each protein
psm.count.table = data.frame(count = as.matrix(data_psm[,2]), row.names =  data_psm$Gene.Symbol)
fit3$count = psm.count.table[rownames(fit3$coefficients),"count"]
fit4 = spectraCounteBayes(fit3)

# n=30 limits the boxplot to show only proteins quantified by <= 30 PSMs.
VarianceBoxplot(fit4,n=100,main="TMT PSM count",xlab="PSM count")
VarianceScatterplot(fit4,main="TMT PSM count")


##Extract outputs from DEqMS (choose one)
#if you are not sure which coef_col refers to the specific contrast,type
head(fit4$coefficients)
DEqMS.results = outputResult(fit4,coef_col = 1)
DEqMS.results = outputResult(fit4,coef_col = 2)
DEqMS.results = outputResult(fit4,coef_col = 3)

# Add Gene names to the data frame
DEqMS.results <- tibble::rownames_to_column(DEqMS.results, "Gene.Symbol")



#### HEATMAP (Figure 4D) ####
group_colors <- list(Group=c(Acute_Control="snow3", Acute_sgEif4g2="lightskyblue", Chronic_Control='black', Chronic_sgEif4g2="blue"))

#You can select the top expressed genes across all samples.
num_rows <- nrow(data.log)
print(num_rows)
select <- order(rowMeans(data.log), decreasing=TRUE)[1:2000]

#Or you can select specific genes
select <- c('Stmn1', 'Cdk1', 'Cdk2', 'Bub1b', 'Mcm3', 'Mcm4', 'Cdca2', 'Ccnb1')

select <- c('Slc7a11', 'Aifm2', 'Alox5', 'Gclc', 'Acsl3', 'Mboat7', 'Prdx2', 'Prdx4', 'Apaf1', 'Pycard')

select <- c('Hspbap1', 'Hspa4l', 'Hspg2')

select <- c('Cd44', 'Il2ra', 'Il12rb2', 'Il18r1', 'Icos',
            'Plcg1', 'Irf4', 'Irf5', 'Tbx21', 'Prdm1', 'Batf3')

select <- c('Pdcd1', 'Tnfrsf9', 'Cd38', 'Cd244', 'Cd69', 'Cd160', 'Snx9', 'Tgfb1', 'Nfatc1',
            'Bhlhe40', 'Ctla4', 'Lag3')

select <- c('Il7r', 'Foxo1', 'Ccr7', 'Cd27')

select <- c('Klra7', 'Klra4',  'Klrk1', 'Klrd1', 'Klra2')

#Make dataframe for group labels
df <- data.frame(Group = rep(c("Acute_Control", "Acute_sgEif4g2", "Chronic_Control", "Chronic_sgEif4g2"), c(3, 3, 3, 3)))
row.names(df) <- colnames(data)

#Convert counts to matrix
data <- as.matrix(data.log)

any(is.na(data))
any(is.nan(data))
any(is.infinite(data))

color_palette <- colorRampPalette(brewer.pal(9, "RdBu"))(100)
reversed_palette <- rev(color_palette)

#The pheatmap function does quite a lot in a single step, including scaling the data by row and clustering both the samples (columns) and genes (rows).
p <- pheatmap(data[select,], name='Expression', scale="row",  cluster_rows=TRUE, show_rownames=FALSE, 
              cluster_cols=FALSE, annotation_col=df, color = reversed_palette,
              labels_col=c('Acute_Control.1', 'Acute_Control.2', 'Acute_Control.3', 'Acute_sgEif4g2.1', 'Acute_sgEif4g2.2', 'Acute_sgEif4g2.3',
                           'Chronic_Control.1', 'Chronic_Control.2', 'Chronic_Control.3', 'Chronic_sgEif4g2.1', 'Chronic_sgEif4g2.2', 'Chronic_sgEif4g2.3'),
              annotation_colors=group_colors, fontsize = 9, fontsize_row = 9, main = 'Proteomics top 2000')
p <- pheatmap(data[select,], name='Expression', scale="row",  
              cluster_rows=FALSE, cluster_cols=FALSE, show_rownames=TRUE, show_row_dend = FALSE, annotation_col=df,
              labels_col=c('Acute_Control.1', 'Acute_Control.2', 'Acute_Control.3', 'Acute_sgEif4g2.1', 'Acute_sgEif4g2.2', 'Acute_sgEif4g2.3',
                           'Chronic_Control.1', 'Chronic_Control.2', 'Chronic_Control.3', 'Chronic_sgEif4g2.1', 'Chronic_sgEif4g2.2', 'Chronic_sgEif4g2.3'),
              annotation_colors=group_colors, color = reversed_palette, fontsize = 10, fontsize_row = 10, 
              cellwidth=10, cellheight = 10, main = 'Proteomics')

dev.off()
p


#Add vertical lines to heatmap
co=column_order(p)
nc = ncol(p)

decorate_heatmap_body(heatmap = "Expression", code = {
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



#### PCA Plot (Figure S4D) ####
library(ggfortify)
#Generate sample annotation table.
cond = factor(c(rep(c("Acute_Control"), 3), rep(c("Acute_sgEif4g2"), 3), 
                rep(c("Chronic_Control"), 3), rep(c("Chronic_sgEif4g2"), 3))) 
sampleTable <- data.frame(
  row.names = colnames(data)[1:12],
  cond = as.factor(cond))

all(colnames(data) == rownames(sampleTable))

#Generate the PCA projections
log_counts_transposed <- t(data.log) #make samples as rows and genes as columns
pc <- prcomp(log_counts_transposed, scale. = FALSE) #scale=TRUE means the data gets normalized!
df <- as.data.frame(pc$x[,1:2]) 
df$PC1 <- as.numeric(df$PC1)
df$PC2 <- as.numeric(df$PC2)

#Merge Sample table with PCA dataframe
df <- transform(merge(sampleTable, df, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)

Group <- factor(df$cond)
SampleName <- factor(row.names(df))

#Generate a customized plot using the ggplot2 package, including shape assignments for our treatment groups and color assignments for the individual samples
#Different shapes for acute and chronic (fill/shape in aes)
gp <- ggplot(df, aes(x = PC1, y = PC2), fill = cond) +
  geom_point(aes(fill=cond), size=5, shape=21, color = 'black') + #Shape refers to what kind of pointshape/outline (21 is outlined circle)
  scale_fill_manual(values = c(Acute_Control = "snow3", Acute_sgEif4g2 = "lightskyblue",
                               Chronic_Control = 'black', Chronic_sgEif4g2 = "blue")) +
  ggtitle(label = as.character('PCA'))  + theme_bw(base_size = 10) +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position.inside = c(0.5,0.5), legend.key.size = unit(0.4, "cm"),
        legend.background = element_rect(fill="white", size=0.5, linetype="solid", colour ="black"))

gp
ggsave('PCA_plot.jpg', plot = gp, width = 4, height = 4)



#### PATHWAY ANALYSIS (Figure S4G) ####
library(tidyverse)
library(org.Mm.eg.db)
library(ReactomePA)
library(msigdbr)
library(ComplexHeatmap)
library(clusterProfiler)

DEqMS.results <- DEqMS.results[order(-DEqMS.results$logFC), ] #Order in decreasing LFC (required for Clusterprofiler)
gene_list <- DEqMS.results$logFC #Make list with just gene name and LFC
names(gene_list) <- DEqMS.results$Gene.Symbol
gene_list<-na.omit(gene_list)

## MSIGDB HALLMARK
gsea_h_df <- data.frame()

m_df <- msigdbr(species = 'Mus musculus', category = 'H') %>% dplyr::select(gs_name, gene_symbol)
gsea_h <- GSEA(gene_list, TERM2GENE = m_df, pAdjustMethod = 'BH', pvalueCutoff = 0.05)
gsea_h_df <- gsea_h@result

gseaplot2(gsea_h, geneSetID = c(7, 8), pvalue_table = FALSE, subplots = 1:2, base_size = 15,
          color = c("palevioletred1", "mediumseagreen")) 
gseaplot2(gsea_h, geneSetID = c(1, 2), pvalue_table = FALSE, subplots = 1:2, base_size = 15,
          color = c("palevioletred1", "mediumseagreen")) 


