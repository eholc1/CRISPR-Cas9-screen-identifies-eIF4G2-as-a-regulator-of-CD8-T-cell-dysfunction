#Single-cell analysis

library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(scCustomize)
library(viridis)

#### Seurat pre-processing ####
#load all sample data
PBS1<- Read10X('./data/PBS_1/')
PBS2<- Read10X('./data/PBS_2/')
PBS3<- Read10X('./data/PBS_3/')

PBS1_SO = CreateSeuratObject(counts = PBS1, project = 'PBS1')
PBS2_SO = CreateSeuratObject(counts = PBS2, project = 'PBS2')
PBS3_SO = CreateSeuratObject(counts = PBS3, project = 'PBS3')

CTRL1<- Read10X('./data/Ctrl_1/')
CTRL2<- Read10X('./data/Ctrl_2/')
CTRL3<- Read10X('./data/Ctrl_3/')

CTRL1_SO = CreateSeuratObject(counts = CTRL1, project = 'CTRL1')
CTRL2_SO = CreateSeuratObject(counts = CTRL2, project = 'CTRL2')
CTRL3_SO = CreateSeuratObject(counts = CTRL3, project = 'CTRL3')

EIF1<- Read10X('./data/Eif_1/')
EIF2<- Read10X('./data/Eif_2/')
EIF3<- Read10X('./data/Eif_3/')

EIF1_SO = CreateSeuratObject(counts = EIF1, project = 'EIF1')
EIF2_SO = CreateSeuratObject(counts = EIF2, project = 'EIF2')
EIF3_SO = CreateSeuratObject(counts = EIF3, project = 'EIF3')

#merge datasets
data <- merge(PBS1_SO, y = c(PBS2_SO, PBS3_SO, CTRL1_SO, CTRL2_SO, CTRL3_SO, EIF1_SO, EIF2_SO, EIF3_SO))

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
#mitochondrial genes
data[['percent.mt']] <- PercentageFeatureSet(data, pattern = '^mt-')

#ribosomal genes
data[['percent.ribo']] <- PercentageFeatureSet(object = data, pattern = '^Rp')

# Visualize QC metrics as a violin plot with no dots shown
VlnPlot(data, features = c('nFeature_RNA'), pt.size = 0, raster = FALSE)
VlnPlot(data, features = c('nCount_RNA'), pt.size = 1, raster = FALSE, y.max = 300000)
VlnPlot(data, features = c('percent.mt'), pt.size = 1, raster = FALSE)
VlnPlot(data, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, pt.size = 1, raster = FALSE)
VlnPlot(data, features = 'nCount_RNA', pt.size = 0, y.max = 6000)

#apply QC measurements
data <- subset(data, subset = nCount_RNA<300000 & nFeature_RNA<8000 & percent.mt<20)


#### Data analysis ####
#normalization and scaling
data <- NormalizeData(data) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
RunPCA(data)

#PCA
DimPlot(data, reduction = 'pca', group.by = 'orig.ident', raster = FALSE)
DimPlot(data, reduction = 'pca', split.by = 'orig.ident', ncol = 3)

#Get dimensions
pct <- data@reductions$pca@stdev / sum(data@reductions$pca@stdev) * 100 
cum <- cumsum(pct) 
co1 <- which(cum > 90 & pct < 5)[1] 
co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > .1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
pcs
#Elbow plot method
ElbowPlot(data)

#Clustering
#We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells
data <- FindNeighbors(data, dims = 1:18)
data <- FindClusters(data, resolution = 1.2)
data <- RunUMAP(data, dims = 1:18)

DimPlot (data, reduction = 'umap', raster = FALSE, label = TRUE, label.size = 5)
DimPlot (data, reduction = 'umap', split.by = 'orig.ident', raster = FALSE, ncol = 3)
DimPlot (data, reduction = 'umap', split.by = 'seurat_clusters', raster = FALSE, ncol = 3)

#Combine triplicates
data <- JoinLayers(data)

saveRDS(data, file = 'data.rds')
data <- readRDS('data.rds')


# Stash cell identity classes(cluster numbers)
data[['replicate.ident']] <- data@meta.data$orig.ident

#Rename idents
data <- SetIdent(data, value = 'orig.ident')
data<- RenameIdents(data, 'PBS1' = 'PBS', 'PBS2' = 'PBS', 'PBS3' = 'PBS', 
                       'CTRL1' = 'CTRL', 'CTRL2' = 'CTRL', 'CTRL3' = 'CTRL',
                       'EIF1' = 'EIF', 'EIF2' = 'EIF', 'EIF3' = 'EIF')

# Stash cell identity genotypes
data[['group.ident']] <- Idents(data)

#set ident to seurat cluters
data <- SetIdent(data, value = 'seurat_clusters')


#Find cluster markers
datamarkers1 <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
datamarkers1  %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
#save as a .csv
write.csv(datamarkers1, file = 'data_markers.csv')

#Rename clusters 
new.cluster.ids <- c('NK', 'NKT', 'CD8/CD4', 'C1q+ Mac', 'Cxcl10+Il18+ mo-Mac', 'CD4', 'Klf4+Tlr7+ mo-Mac', 
                     'B cell', 'CD8', 'CD8', 'CD8', 'Arg1+ mo-Mac', 'CD8',
                     'CD8', 'CD4', 'C1q+Trem2+ Mac', 'Melanocyte', 'Neutrophil', 'mregDC', 'pDC', 'cDC1', 
                     'CD8', 'Ifng+Ccl5+ mo-Mac', 'cDC2', 'CD4', 'Melanocyte', 'NK', 'CD8', 'RBC', 'B cell',
                     'CD8', 'Mast cell', 'B cell')

names(new.cluster.ids) <- levels(data)
data_id <- RenameIdents(data, new.cluster.ids)

levels(data_id) <- c('NK', 'NKT', 'CD8', 'CD4', 'CD8/CD4', 'B cell',
                    'C1q+ Mac', 'C1q+Trem2+ Mac', 'Ifng+Ccl5+ mo-Mac', 'Cxcl10+Il18+ mo-Mac', 'Klf4+Tlr7+ mo-Mac', 'Arg1+ mo-Mac', 
                    'mregDC', 'cDC1', 'cDC2', 'pDC', 'Neutrophil', 'Mast cell', 'Monocyte/Epi', 'Melanocyte', 'RBC')
data_id[['cell.ident']] <- Idents(data_id)

#saveRDS
saveRDS(data_id, file = 'data_id.rds')
data_id <- readRDS('data_id.rds')


#### UMAP, Frequency plot, & Dotplot (Figure S2C-E) ####
#UMAP
mypalette<-DiscretePalette_scCustomize(num_colors=20, palette='varibow', shuffle=FALSE)
DimPlot(data_id, label = TRUE, label.size = 4, raster = FALSE, cols=mypalette)
DimPlot(data_id, label = TRUE, label.size = 4, split.by = 'group.ident', raster = FALSE)

#Frequency plot
mypalette<-DiscretePalette_scCustomize(num_colors=20, palette='varibow', shuffle=FALSE)

x <- data.frame(table(data_id@meta.data$group.ident, data_id@active.ident))
x$Var2 <- factor(x$Var2)
x <- reshape2::dcast(x, Var1~Var2) 
rownames(x) <- x$Var1
x <- x[,-1]
x <- x/rowSums(x)
x$Var1 <- rownames(x)
x <- x %>% as.data.frame()
x <- reshape2::melt(x, ids.vars = colnames(x))
x$Var1 <- factor(x$Var1, levels = c('PBS', 'CTRL', 'EIF'))
x <- na.omit(x)
print(x %>% ggplot(aes(x = Var1, y = value, fill = variable)) + geom_bar(position = 'stack', 
                                                                         stat = 'identity') + theme_classic() + labs(x = 'Treatment', y = 'Percentage') + 
        theme(legend.title = element_blank())) + scale_fill_manual(values = mypalette)

#Dotplot defining clusters
DotPlot(data_id, features = c('Ncr1', 'Klrb1c', 'Klra6', 'Trdc', 
                              'Cd3e', 'Cd8a', 'Cd4', 
                              'Cd19', 'Cd22', 
                              'Cd68',  'Ccr2',
                              'C1qa', 'C1qb',
                              'Trem2', 'Mrc1', 'Ccl8',
                              'Ifng', 'Ccl5',
                              'Ccl2', 'Cxcl10', 'Il18', 'Tnf',
                              'C3', 'Lyz1', 'Tlr7', 'Klf4',
                              'Arg1', 'Nos2', 'Cd274',
                              'Flt3', #pan DC
                              'Etv3', 'Ccr7', 'Fscn1', 'Cd80', 'Il12b', 'Mreg',  #migratory DC?
                              'Irf8', 'Xcr1', 'Itgae', #cDC1
                              'Zbtb46', 'Irf4', 'Cd209a', #cDC2
                              'Tlr9', 'Bst2', 'Siglech',  #pDC
                              'S100a9', 'Mmp9', 'Cxcr2',
                              'Mcpt8', 'Il6', 'Il4', 'Ccl3',
                              'Pmel', 'Mlana'
                              )) + scale_colour_gradientn(colours = rev(brewer.pal(11, "RdYlBu"))) + 
                                   theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank()) + RotatedAxis()


#### CD8T subsetting and filtering ####
subset_CD8 = subset(data_id, idents = c('CD8','CD8/CD4'))
subset_CD8 <- SCTransform(object = subset_CD8, verbose = FALSE)
subset_CD8 <- RunPCA(object = subset_CD8, verbose = FALSE)
DimPlot(subset_CD8, reduction = 'pca', split.by = 'orig.ident')
DimPlot(subset_CD8, reduction = 'pca')

pct <- subset_CD8@reductions$pca@stdev / sum(subset_CD8@reductions$pca@stdev) * 100 
cum <- cumsum(pct) 
co1 <- which(cum > 90 & pct < 5)[1] 
co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > .1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
pcs

ElbowPlot(subset_CD8)
subset_CD8 <- FindNeighbors(object = subset_CD8, dims = 1:14, verbose = FALSE)
subset_CD8 <- FindClusters(subset_CD8, resolution = 1.5)
subset_CD8 <- RunUMAP(object = subset_CD8, dims = 1:14, verbose = FALSE)
mypalette<-DiscretePalette(16, palette='glasbey', shuffle=FALSE)
DimPlot(subset_CD8, reduction = 'umap', label = TRUE)
DimPlot(subset_CD8, reduction = 'umap', cols=mypalette, label = TRUE)
DimPlot(subset_CD8, reduction = 'umap', split.by = 'orig.ident', label = FALSE)

#markers for new Tcell clusters
CD8_markers <- FindAllMarkers(subset_CD8, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(CD8_markers, file = 'subset_CD8.csv')

#Rename clusters
CD8.clusters <- c('Cycling', 'Naive', 'Fgl2+Ccl5+','CD4', 'Cxcr3+ TRM', 'Fgl2+Ccl5+', 'Cycling', 
                  'Terminal Effector', 'Gzma+ TEM', 'Trbv29+', 'Mac', 'CD4', 'ISAG', 'Gzm+', 
                  'Naive', 'Il2ra+Xcl1+', 'Mito Unident', 'Naive', 'Cycling', 'Naive', 
                  'Naive', 'Cycling', 'HSP+', 'Fgl2+Ccl5+', 'Naive',
                  'Mac', 'CD4')

CD8.ids <- SetIdent(CD8.ids, value = 'SCT_snn_res.1.5')
names(CD8.clusters) <- levels(CD8.ids)
CD8.ids <- RenameIdents(CD8.ids, CD8.clusters)

levels(CD8.ids) <- c('Naive', 'Cxcr3+ TRM', 'Gzma+ TEM', 
                     'Il2ra+Xcl1+', 'Gzm+',
                     'Terminal Effector', 'HSP+', 'Trbv29+', 'Fgl2+Ccl5+', 
                     'ISAG', 'Cycling', 'Mito Unident',
                     'Mac', 'CD4')

DimPlot(CD8.ids, label=FALSE)
DimPlot(CD8.ids, reduction = 'umap', split.by = 'group.ident', label = FALSE, cols=mypalette)


#Save RDS
saveRDS(subset_CD8, file = 'subset_CD8.rds')
subset_CD8 <- readRDS("subset_CD8.rds")    
saveRDS(CD8.ids, file = 'CD8.ids.rds')
CD8.ids <- readRDS("CD8.ids.rds")    



#Filter out residual CD4 and Macrophage from CD8 T cells
CD8_filter = subset(CD8.ids, idents = c('Naive', 'Cxcr3+ TRM', 'Gzma+ TEM', 
                                        'Il2ra+Xcl1+', 'Gzm+',
                                        'Terminal Effector', 'HSP+', 'Trbv29+', 'Fgl2+Ccl5+', 
                                        'ISAG', 'Cycling', 'Mito Unident'))


levels(CD8_filter) <- c('Naive', 'Cxcr3+ TRM', 'Gzma+ TEM', 
                        'Il2ra+Xcl1+', 'Gzm+',
                        'Terminal Effector', 'HSP+', 'Trbv29+', 'Fgl2+Ccl5+', 
                        'ISAG', 'Cycling', 'Mito Unident')

#Rename clusters
CD8.clusters <- c('Naive', 'Cxcr3+ TRM', 'Gzma+ TEM', 
                  'Il2ra+Xcl1+', 'Gzm+',
                  'Terminal Effector', 'HSP+', 'Trbv29+', 'Fgl2+Ccl5+', 
                  'ISAG', 'Cycling', 'Mito Unident')

names(CD8.clusters) <- levels(CD8_filter)
CD8_filter <- RenameIdents(CD8_filter, CD8.clusters)

#markers for new clusters
tcell_markers <- FindAllMarkers(CD8_filter, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(tcell_markers, file = 'CD8_filtered.csv')

saveRDS(CD8_filter, file = 'CD8_filter.rds')
CD8_filter <- readRDS("CD8_filter.rds")


#### CD8T UMAP & Frequency plot (Figure 2H-I) ####
mypalette<-DiscretePalette_scCustomize(num_colors=12, palette='varibow', shuffle=FALSE)
DimPlot(CD8_filter, label=FALSE, cols=mypalette)
DimPlot(CD8_filter, reduction = 'umap', split.by = 'group.ident', label = FALSE, cols=mypalette)

#Frequency plot
x <- data.frame(table(CD8_filter@meta.data$group.ident, CD8_filter@active.ident))
x$Var2 <- factor(x$Var2)
x <- reshape2::dcast(x, Var1~Var2) 
rownames(x) <- x$Var1
x <- x[,-1]
x <- x/rowSums(x)
x$Var1 <- rownames(x)
x <- x %>% as.data.frame()
x <- reshape2::melt(x, ids.vars = colnames(x))
x$Var1 <- factor(x$Var1, levels = c('PBS', 'CTRL', 'EIF'))
x <- na.omit(x)
print(x %>% ggplot(aes(x = Var1, y = value, fill = variable)) + geom_bar(position = 'stack', 
                                                                         stat = 'identity') + theme_classic() + labs(x = 'Treatment', y = 'Percentage') + 
        theme(legend.title = element_blank())) + scale_fill_manual(values = mypalette)

#### CD8T Dotplot (Figure S2F) ####
DotPlot(CD8_filter, features = c('Satb1', 'Bach2', 'Sell', 'Ccr7', 'Lef1', 'Tcf7', 'S1pr1',
                                 'Il7r', 'Cd44', 'Abcb1a', 'Cxcr3', 
                                 'Gzma', 'Ccr2',  'Itga1',
                                 'Xcl1', 'Ccl1', 'Cd160', 'Il2ra',
                                 'Gzmc', 'Gzmd', 'Gzme', 'Gzmf', 'Ccl9', 'Cd244a',
                                 'Ccl3', 'Ccl4', 'Ifng', 'Tnf', 'Prf1', 'Havcr2', 'Rgs16', 'Lag3',
                                 'Tox', 'Tigit',  'Rgs1', 'Ctla4', 'Prdm1',
                                 'Hspa1a', 'Hspa1b', 'Dnajb1',
                                 'Klri2', 'Trbv29','Klrd1',
                                 'Fgl2', 'Ccl5', 'Ccr5',
                                 'Ifit1', 'Isg15', 'Cxcl10', 
                                 'Mki67', 'Cdk1', 'Ska1'
)) +  theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank()) + scale_colour_gradientn(colours = rev(brewer.pal(11, "RdYlBu"))) + RotatedAxis()


#### CD8T Violin Plots (Figures 2K-N and S2G-I) ####
tcell_sig <- read.csv("T_cell_signatures.csv", header=TRUE)
tcell_sig <- tcell_sig %>% dplyr::select(GO_tcell_prolif)  #Select which signature
tcell_sig <- tcell_sig[!duplicated(tcell_sig), ] # all duplicate genes removed
tcell_sig <- data.frame(tcell_sig)
tcell_sig <- tcell_sig[!(tcell_sig ==""), ]
tcell_sig <- data.frame(tcell_sig)
tcell_sig <- tcell_sig %>% dplyr::rename(Gene = tcell_sig)
tcell_sig$Gene <- trimws(tcell_sig$Gene)

## Convert human to mouse gene names if applicable
  mouse_human_genes <- read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
  
  # separate human and mouse 
  mouse <- split.data.frame(mouse_human_genes,mouse_human_genes$Common.Organism.Name)[[2]]
  human <- split.data.frame(mouse_human_genes,mouse_human_genes$Common.Organism.Name)[[1]]
  
  # remove some columns
  mouse <- mouse[,c(1,4)]
  human <- human[,c(1,4)]
  
  # merge the 2 dataset  (note that the human list is longer than the mouse one)
  mh_data <- merge.data.frame(mouse,human,by = "DB.Class.Key",all.y = TRUE) 
  
  #Rename column to whatever the column is you want to merge with
  names(mh_data)[3] <- paste("Gene")
  
  #Join files based on gene set list
  tcell_sig <- dplyr :: left_join(tcell_sig, mh_data)
  tcell_sig <- dplyr :: select(tcell_sig, Symbol.x)
  tcell_sig <- na.omit(tcell_sig)

#Add module score
CD8_filter <- AddModuleScore(CD8_filter, features = tcell_sig, ctrl = 100,  name = 'GO_tcell_prolif')

##Violin Plot
group_colors <- c("PBS" = "grey94", "CTRL" = "grey60", "EIF" = "blue")
#select comparison
my_comparisons <- list( c("CTRL", "PBS"))
my_comparisons <- list( c("EIF", "CTRL"))
my_comparisons <- list(c("EIF", "CTRL"), c("CTRL", "PBS"))

VlnPlot(subset(CD8_filter), features = c('Pancan_ISG1'), 
  ncol = 1, group.by = 'group.ident', pt.size = 0, cols = group_colors) + 
  NoLegend() + theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.title.x = element_blank()) + 
  scale_y_continuous(expand = expansion(mult = c(0., 0.1))) + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +  
  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)

ggsave(width = 4, height = 4, dpi = 300, "plot.png")


# Extract module scores
module_scores <- data.frame(
  group = CD8_filter$group.ident,
  score = CD8_filter$Pancan_ISG1)

filtered_scores <- module_scores %>%
  filter(group %in% c("CTRL", "EIF"))

# Wilcoxon Test for two groups
wilcox_test_result <- wilcox.test(score ~ group, data = filtered_scores)
print(wilcox_test_result)


#### GSEA Analysis: MSIGDB Hallmark (Figure 2J) ####
library(tidyverse)
library(Seurat)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)
library(msigdbr)
library(ComplexHeatmap)

ob <- CD8_filter
idents_1 <- c('EIF')
idents_2 <- c('CTRL')
gsea_h_df <- data.frame()

for (i in 1:length(idents_1)) {
  
  cat('testing = ', idents_1[i], ' vs. ', idents_2[i], '\n')
  
  fc <- FoldChange(ob, ident.1 = idents_1[i], ident.2 = idents_2[i], group.by = 'group.ident')
  
  bitr_df <- bitr(rownames(fc), fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
  
  fc$ENTREZID <- bitr_df$ENTREZID[match(rownames(fc),bitr_df$SYMBOL)]
  fc <- fc[!is.na(fc$ENTREZID),]
  
  gene_list <- fc %>% filter(abs(avg_log2FC) > 0) %>% arrange(-avg_log2FC) %>% dplyr::select(ENTREZID,avg_log2FC) %>% tibble::deframe()
  
  
  ## MSIGDB
  ### H
  m_df <- msigdbr(species = 'Mus musculus', category = 'H') %>% dplyr::select(gs_name, entrez_gene)
  gsea_h <- GSEA(gene_list, TERM2GENE = m_df, pAdjustMethod = 'BH', pvalueCutoff = 0.05)
  gsea_h@result$comparison <- paste0(idents_1[i],'_',idents_2[i])
  gsea_h@result$core_enrichment <- sapply(gsea_h@result$core_enrichment, function(i) {
    y <- strsplit(i,'/')[[1]]
    y <- bitr_df$SYMBOL[match(y,bitr_df$ENTREZID)]
    paste(y, collapse = '/')
  })
  rownames(gsea_h@result) <- NULL
  gsea_h_df <- rbind(gsea_h_df, gsea_h@result)
  gc()
  
}

write_csv(gsea_h_df, file = 'cd8_gsea_msigdb_hallmark.csv')


#GSEA Plot
#You can select which gene set to view. The first gene set is 1, second gene set is 2, etc. Default: 1  
gseaplot(gsea_h, by = 'all', title = gsea_h$Description[1], geneSetID = 1)
dev.off()

