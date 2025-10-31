#Ribosome profiling analysis

library(ggplot2)
library(ggrepel)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(ComplexHeatmap)
library(matrixStats)
library(limma)
library(tibble)

rm(list = ls())

#### Data Import ####
#Folder location: RNAseq -> Salmon.STD_TS
data = read.csv("gene_agg_cpm_log_matrix.csv", header=T, row.names=1) #RNA data
#Folder location: Riboseq -> Salmon.STD_TS.STD_TL
data = read.csv("gene_agg_cpm_log_matrix.csv", header=T, row.names=1) #Ribo data

#Reorder columns ONLY FOR RIBO DATA
data <- data[, c(4:12, 1:3)]

#Remove NA values
data <- data %>% replace(is.na(.), 0)

#Add gene symbol to dataframe
library(org.Mm.eg.db)
library(clusterProfiler)
bitr_df <- bitr(rownames(data), fromType = 'ENSEMBL', toType = 'SYMBOL', OrgDb = org.Mm.eg.db)
data$SYMBOL <- bitr_df$SYMBOL[match(rownames(data),bitr_df$ENSEMBL)]
data <- data[!is.na(data$SYMBOL),]
data <- rownames_to_column(data, var = "ENSEMBL") #Make row names into column

#Check for duplicate genes
has_duplicates <- any(duplicated(data$SYMBOL))
print(has_duplicates)

# Deal with duplicated gene symbols by combining gene symbol with ENSEMBL ID if non-unique
data$SYMBOL <- ifelse(
  duplicated(data$SYMBOL),
  paste(data$SYMBOL, data$ENSEMBL, sep="_"),
  data$SYMBOL)

rownames(data) <- data$SYMBOL #Make gene the row names
data <- data[, !(names(data) %in% c("ENSEMBL", "SYMBOL"))] #Make dataframe with just counts



#### Riboseq PCA Plot (Figure S4I) ####
library(ggfortify)
#Generate sample annotation table.
cond = factor(c(rep(c("Acute_Control"), 3), rep(c("Acute_sgEif4g2"), 3), 
                rep(c("Chronic_Control"), 3), rep(c("Chronic_sgEif4g2"), 3))) 
sampleTable <- data.frame(
  row.names = colnames(data)[1:12],
  cond = as.factor(cond))

all(colnames(data) == rownames(sampleTable))

#Generate the PCA projections
data.transposed <- t(data) #make samples as rows and genes as columns
pc <- prcomp(data.transposed, scale. = FALSE) #scale=TRUE means the data gets normalized!
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
  ggtitle(label = as.character('RNAseq'))  + theme_bw(base_size = 10) +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position.inside = c(0.5,0.5), legend.key.size = unit(0.4, "cm"),
        legend.background = element_rect(fill="white", size=0.5, linetype="solid", colour ="black"))

gp
ggsave('PCA_plot.jpg', plot = gp, width = 4, height = 4)



#### Riboseq Heatmap (Figure S4J) ####
#Heatmap with top 500 variant or expressed genes, rlog normalized data
group_colors <- list(Group=c(Acute_Control="snow3", Acute_sgEif4g2="lightskyblue", Chronic_Control='black', Chronic_sgEif4g2="blue"))

#Select the top expressed genes across all samples 
num_rows <- nrow(data)
print(num_rows)
select <- order(rowMeans(data), decreasing=TRUE)[1:2000]

#Make dataframe for group labels
df <- data.frame(Group = rep(c("Acute_Control", "Acute_sgEif4g2", "Chronic_Control", "Chronic_sgEif4g2"), c(3, 3, 3, 3)))
row.names(df) <- colnames(data)

#Convert counts to matrix
data <- as.matrix(data)

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
              annotation_colors=group_colors, fontsize = 9, fontsize_row = 9, main = 'Riboseq top 2000')

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



#### Riboseq Volcano Plot (Figures 4E-F) ####
rm(list = ls())

#Select Riboseq differential expression (DE) file for comparison
data = read.csv("CD8T_Eif4g2_KO_and_ACU_AGST__CD8T_and_ACU_AGST_Ribo.csv", header=T, row.names=1) #Acute KO vs. Acute control
data = read.csv("CD8T_and_CHR_AGST__CD8T_and_ACU_AGST_Ribo.csv", header=T, row.names=1) #Chronic control vs. Acute control
data = read.csv("CD8T_Eif4g2_KO_and_CHR_AGST__CD8T_and_CHR_AGST_Ribo.csv", header=T, row.names=1) #Chronic KO vs. Chronic control
#Change ENSEMBL to SYMBOL (code located in data import section above)

Comparison <- "Acute sgEif4g2 vs. Acute Control"
Comparison <- "Chronic Control vs. Acute Control"
Comparison <- "Chronic sgEif4g2 vs. Chronic Control"

fc <- 1.5
pval <- 0.05

#Subset data to label the datapoints (genes) that pass our thresholds.
data$dot <- rep(3, nrow(data))
data$dot[which(data$FDR <= pval & data$logFC < 0 & abs(data$logFC) >= log2(fc))] = 2
data$dot[which(data$FDR <= pval & data$logFC > 0 & abs(data$logFC) >= log2(fc))] = 1
# 1 = UP, 2= DOWN

#Write csv to view all DE genes
write.csv(data, file="Ribo.AE-AC.csv")
write.csv(data, file="Ribo.CC-AC.csv")
write.csv(data, file="Ribo.CE-CC.csv")


#Volcano Plot
#Acute sgEif4g2 vs. Acute Control
select <- c('Eif4g2', 
            'Cd69', 'Zap70', 'Cd3g', 'Jund', 'Jun', 'Akt2', 'Stat6', 'Nfatc1', 'Nfatc3', 'Lat', 'Dusp6', 'Dusp1', 'Dusp10', 'Dusp2', 'Fyn', #Activation
            'Gzmc', 'Gzmd', 'Gzmf', 'Tnf', 'Prf1', 'Il16', 'Klrg1', 'Myc', 'Max', 'Maf', 'Cxcl10', 'Ccl1', #Effector
            'Fli1', 'Cd101', 'Cd160', #Exhaustion
            'Atf6', 'Atf4', 'Aim2', 'Card11', #Stress
            'Slc7a11', 'Alox5') #Ferroptosis
gene_groups <- list(
  Activation = c('Cd69','Zap70','Cd3g','Jund','Jun','Akt2','Stat6',
                 'Nfatc1','Nfatc3','Lat','Dusp6','Dusp1','Dusp10','Dusp2','Fyn'),
  Effector   = c('Gzmc','Gzmd','Gzmf','Tnf','Prf1','Il16','Klrg1','Myc','Max','Maf','Cxcl10','Ccl1'),
  Exhaustion = c('Fli1','Cd101','Cd160'),
  Stress     = c('Atf6','Atf4','Aim2','Card11','Slc7a11','Alox5'),
  Etc = c('Eif4g2')
)

           
#Chronic sgEif4g2 vs. Chronic Control
select <- c('Eif4g2', #Translation
            'Cd3g', 'Jun', 'Akt2', 'Sos1', 'Lat', 'Itk', 'Fyn', 'Inpp1', 'Nr4a2', 'Dusp1', 'Dusp6', 'Dusp10', #TCR
            'Gzmc', 'Gzmg', 'Gzmf', 'Prf1', 'Ccl5', 'Tnf', 'Il16', 'Id2', 'Max', 'Klrg1', #Effector
            'Pdcd1',  'Fli1', 'Snx9', 'Cd244a', 'Prdm1', #Exhaustion
            'Atf4', 'Atf5', 'Aim2', 'Pycard', 'Casp3', 'Nlrp3', #Death
            'Slc7a11', 'Alox5') #Ferroptosis
gene_groups <- list(
  Activation = c('Cd3g', 'Jun', 'Akt2', 'Sos1', 'Lat', 'Itk', 'Fyn', 'Inpp1', 'Nr4a2', 'Dusp1', 'Dusp6', 'Dusp10'),
  Effector   = c('Gzmc', 'Gzmg', 'Gzmf', 'Prf1', 'Ccl5', 'Tnf', 'Il16', 'Id2', 'Max', 'Klrg1'),
  Exhaustion = c('Pdcd1',  'Fli1', 'Snx9', 'Cd244a', 'Prdm1'),
  Stress     = c('Atf4', 'Atf5', 'Aim2', 'Pycard', 'Casp3', 'Nlrp3', 'Slc7a11', 'Alox5'),
  Etc = c('Eif4g2')
)


#Create a lookup table
group_df <- stack(gene_groups)
colnames(group_df) <- c("gene_name","group")

#Merge with your volcano data
#Now `data$group` contains the functional group for select genes, NA for others
#You can map color or shape to this new variable
data <- merge(data, group_df, by="gene_name", all.x=TRUE)

#Split data by direction for labeling
label_up <- subset(data, gene_name %in% select & logFC > 0)
label_down <- subset(data, gene_name %in% select & logFC < 0)

data$genelabel <- factor(data$gene_name, levels = select)

#Count the number of significant up and down genes, assign value for legend
data$dot <- factor(data$dot,levels = c(1,2,3), labels = c(paste0('Up: ', sum(data$dot == 1)),paste0('Down: ', sum(data$dot == 2)),'NS'))

#Custom palette for labels
group_colors <- c(
  Activation = "maroon1",
  Effector   = "springgreen4",
  Exhaustion = "blue",
  Stress     = "#FF7F00",
  Etc = "black"
)

#Generating volcano plot
rm(p) 
p <- ggplot(data, aes(x = logFC, y = -log10(FDR))) + 
  geom_point(aes(color = data$dot), size = .5) + theme_classic() + xlab('Log2 fold-change') + ylab('-log10(p.adj)') + 
  scale_color_manual(name = '', values=c('#B31B21', '#1465AC', 'darkgray'), guide = 'none') + 
  geom_vline(xintercept = c(0, -log2(fc), log2(fc)), linetype = c(1, 2, 2), color = c('black', 'black', 'black')) + geom_hline(yintercept = -log10(pval), linetype = 2, color = 'black') +
  coord_flip() +   ggtitle(as.character(Comparison)) +
  ggnewscale::new_scale_color() +
  geom_label_repel(data = label_up, 
                   aes(x = logFC, y = -log10(FDR), label = gene_name, color = group), 
                   fill = "white", segment.color = "black",
                   nudge_x = 5, nudge_y = 1,
                   size = 3, label.size = 0.2, label.padding = 0.2,    
                   force = 1.5, max.overlaps = 50, min.segment.length = 0, segment.alpha = 1) +
  geom_label_repel(data = label_down, 
                   aes(x = logFC, y = -log10(FDR), label = gene_name, color = group), 
                   fill = "white", segment.color = "black",
                   nudge_x = -5, nudge_y = 1,
                   size = 4, label.size = 0.2, label.padding = 0.2,    
                   force = 1.5, max.overlaps = 50, min.segment.length = 0, segment.alpha = 1) +
  scale_color_manual(values = group_colors, name = "Group") +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 3)))  
  

print(p)
dev.off()



#### dTE uORF analysis (Figure 4G) ####
#Folder location: Multi-seq -> DE -> STD_TS.STD_TL -> gene_agg -> dTE
#Choose comparison
data = read.csv("CD8T_Eif4g2_KO_and_ACU_AGST__CD8T_and_ACU_AGST_dTE.csv", header=T) #Acute KO vs. Acute control
data = read.csv("CD8T_and_CHR_AGST__CD8T_and_ACU_AGST_dTE.csv", header=T, row.names=1) #Chronic control vs. Acute control
data = read.csv("CD8T_Eif4g2_KO_and_CHR_AGST__CD8T_and_CHR_AGST_dTE.csv", header=T) #Chronic KO vs. Chronic control

fc <- 1.5
pval <- 0.05

#Subset data to label the datapoints (genes) that pass our thresholds.
data$dot <- rep(3, nrow(data))
data$dot[which(data$FDR <= pval & data$logFC < 0 & abs(data$logFC) >= log2(fc))] = 2
data$dot[which(data$FDR <= pval & data$logFC > 0 & abs(data$logFC) >= log2(fc))] = 1
# 1 = UP, 2= DOWN

#Read in uORF table
df = read.csv("uORF_genes.csv", header=T) 

df$gene_name <- sapply(strsplit(df$Mouse_ORF_ID, ":"), "[", 1)
df <- df[!duplicated(df$gene_name), ] # duplicate genes removed

#Merge tables
merged_df <- left_join(data, df, by = "gene_name")

#Count the number of significant up and down genes, assign value for legend
merged_df$dot <- factor(merged_df$dot,levels = c(1,2,3), labels = c(paste0('Up: ', sum(merged_df$dot == 1)),paste0('Down: ', sum(merged_df$dot == 2)),'NS'))

write.csv(merged_df, file="uORF.AE-AC.csv")
write.csv(merged_df, file="uORF.CE-CC.csv")

