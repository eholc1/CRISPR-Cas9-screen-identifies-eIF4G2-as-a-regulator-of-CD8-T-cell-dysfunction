#Untargeted Lipidomic Analysis
#Targeted Oxidized Phospholipid Analysis

library(dplyr)
library(RColorBrewer)
library(ggrepel)
library(biomaRt)
library(pheatmap)
library(colorspace)
library(ggplot2)

rm(list = ls())

#### DATA PROCESSING ####
data = read.csv("Targeted_data.csv", header=T, row.names=1)
data = read.csv("Untargeted_data.csv", header=T)
data = data[, c(1,10:21)]
data = read.csv("Metaboanalyst_Untargeted.csv", header=T)

data <- data %>% replace(is.na(.), 0)

## Remove duplicate row names for Untargeted file ##
data <- data %>%
  mutate_at(vars(Acute_Control_1, Acute_Control_2, Acute_Control_3,
                 Acute_KO_1, Acute_KO_2, Acute_KO_3,
                 Chronic_Control_1, Chronic_Control_2, Chronic_Control_3,
                 Chronic_KO_1, Chronic_KO_2, Chronic_KO_3), as.numeric)

#Calculate sum of each row, and keep row with highest sum
data$Total <- rowSums(data[, c("Acute_Control_1", "Acute_Control_2", "Acute_Control_3",
                               "Acute_KO_1", "Acute_KO_2", "Acute_KO_3",
                               "Chronic_Control_1", "Chronic_Control_2", "Chronic_Control_3",
                               "Chronic_KO_1", "Chronic_KO_2", "Chronic_KO_3")])
data <- data %>%
  group_by(Acyl.Chain.Level) %>%
  filter(Total == max(Total)) %>%
  ungroup()

duplicates <- data$Acyl.Chain.Level[duplicated(data$Acyl.Chain.Level)] #Check for duplicates
print(duplicates) #Should be 0

data <- as.data.frame(data)
row.names(data) = data$Acyl.Chain.Level
data$Acyl.Chain.Level <- NULL
data$Total <- NULL
data <- data %>% replace(is.na(.), 0)


#Remove rows with all 0's
data <- data[rowSums(data != 0, na.rm = TRUE) > 0, ]

#Normalize to cell number input
data$Acute_Control_1 <- data$Acute_Control_1 / 7.5
data$Acute_Control_2 <- data$Acute_Control_2 / 6.5
data$Acute_Control_3 <- data$Acute_Control_3 / 4.5
data$Acute_KO_1 <- data$Acute_KO_1 / 10
data$Acute_KO_2 <- data$Acute_KO_2 / 10
data$Acute_KO_3 <- data$Acute_KO_3 / 8
data$Chronic_Control_1 <- data$Chronic_Control_1 / 2.7
data$Chronic_Control_2 <- data$Chronic_Control_2 / 2.2
data$Chronic_Control_3 <- data$Chronic_Control_3 / 2.5
data$Chronic_KO_1 <- data$Chronic_KO_1 / 7.86
data$Chronic_KO_2 <- data$Chronic_KO_2 / 6.2
data$Chronic_KO_3 <- data$Chronic_KO_3 / 5.5

data <- as.data.frame(data)


#### PIE CHART (Figure S5A) ####
library(wesanderson)
library(Seurat)
mypalette<-wes_palette("Zissou1", 16, type = "continuous")

#By category
ggplot(data,aes(x= factor(1), fill= factor(data[,"Category"])))+
  geom_bar(stat="count", width=1, color="black") + 
  geom_text(aes(x=1.67, label=scales::percent(..count.. / sum(..count..))), stat = "count", position = position_stack(vjust = 0.5), size=3) +
  coord_polar(theta='y', start = 0, direction = -1) +
  theme_void() + 
  labs(fill = "Category") +
  scale_fill_manual(values = mypalette)

#By Main Class
ggplot(data,aes(x= factor(1), fill= factor(data[,"Main.Class"])))+
  geom_bar(stat="count", width=1, color="black") + 
  geom_text(aes(x=1.67, label=scales::percent(..count.. / sum(..count..))), stat = "count", position = position_stack(vjust = 0.5), size=3) +
  coord_polar(theta='y', start = 0, direction = -1) +
  theme_void() + 
  labs(fill = "Main.Class") +
  scale_fill_manual(values = mypalette)


#### LOAD METABOANALYST ####
devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = FALSE)

metanr_packages <- function(){
  metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest","RBGL","edgeR","fgsea", 'devtools', 'crmn', "httr","qs")
  list_installed <- installed.packages()
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  if(length(new_pkgs)!=0){if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install(new_pkgs)
    print(c(new_pkgs, " packages added..."))
  }
  
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}

library(pacman)
pacman::p_load("iheatmapr", "ellipse", "vegan", "Cairo", "Biobase", "SSPA", "jsonlite", "pkgbuild", "curl","impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", "SSPA", "sva", "limma", "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest","RBGL","edgeR","fgsea", 'devtools', 'crmn', "httr","qs")
library(MetaboAnalystR)



#### METABOANALYST DATA PREP ####
#Prepare Untargeted data for input (output file is supplied)
data = read.csv("Untargeted_data.csv", header=T)
df <- t(data)
df <- as.data.frame(df)
df$Group <- c('Acute_Control', 'Acute_Control', 'Acute_Control', 'Acute_sgEif4g2', 'Acute_sgEif4g2', 'Acute_sgEif4g2',
              'Chronic_Control', 'Chronic_Control', 'Chronic_Control', "Chronic_sgEif4g2", "Chronic_sgEif4g2", "Chronic_sgEif4g2")
df <- df %>%
  tibble::rownames_to_column("ID")
df <- df %>%
  select(ID, Group, everything())

write.csv(df, file="Metaboanalyst_Untargted.csv")
#Now manually delete all rows except for those with count data


#### METABOANALYST DATA INPUT ####
rm(list = ls())
mSet<-InitDataObjects("conc", "stat", FALSE);
mSet<-Read.TextData(mSet, "Metaboanalyst_Targeted.csv", "rowu", "disc");
mSet<-Read.TextData(mSet, "Metaboanalyst_Untargeted.csv", "rowu", "disc");
mSet<-SanityCheckData(mSet)

#Normalize
mSet<-ReplaceMin(mSet)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "LogNorm", "AutoNorm", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_untargeted_", "png", 300, width=NA)
mSet<-PlotSampleNormSummary(mSet, "samplenorm_untargeted", "png", 300, width=NA)

#### METABOANALYST PCA (Figure S5B & S5D) ####
mSet<-PCA.Anal(mSet)
mSet<-GetGroupNames(mSet, "null")
colVec<-c("##919292","##50bcff","##0b0b0b","##0616f3")
shapeVec<-c(0,0,0,0)
mSet<-UpdateGraphSettings(mSet, colVec, shapeVec)

mSet<-PlotPCA2DScore(mSet, "PCA_untargeted", "png", dpi=300, width=7, 1,2, reg=0.95, show=0, grey.scale=0, "na") #With 95% confidence regions
mSet<-PlotPCA2DScore(mSet, "PCA_targeted", "png", dpi=300, width=7, 1,2, reg=0.95, show=0, grey.scale=0, "na") #With 95% confidence regions

#### METABOANALYST HEAT MAP (Figure 6C, S5C, S5E) ####
library(ComplexHeatmap)
#Color the groups
mSet<-GetGroupNames(mSet, "null")
colVec<-c("##919292","##50bcff","##0b0b0b","##0616f3")
shapeVec<-c(0,0,0,0)
mSet<-UpdateGraphSettings(mSet, colVec, shapeVec)

#Can select certain lipids and filter data
data_matrix <- mSet$dataSet$norm
data_matrix <- as.data.frame(data_matrix)

#Glycerophospholipids
selected_lipids <- c('LPC(20:4)', 'LPC(O-16:1)', 'PC(10:0_18:1)', 'PC(14:0_16:1)', 'PC(16:1_16:1)', 'PC(16:1_18:1)', 'PC(26:0)', 'PC(27:0)',
                     'PC(28:1)', 'PC(30:1)', 'PC(30:2)', 'PC(32:2)', 'PC(44:1)', 'LPC(24:1)', 'LPC(O-18:1)', 'PC(18:1_20:4)', 'PC(18:1_24:0)', 
                     'PC(O-31:0)', 'PC(P-16:0_16:0)', 'PC(15:0_20:3)', 'PC(18:1_18:1)', 'PC(26:1)', 'PC(34:3)', 'PC(36:2)', 'PC(O-16:0_14:0)', 
                     'PC(O-16:0_14:1)', 'PC(O-16:1_20:4)', 'PC(P-14:0_18:1)')
selected_lipids <- c('LPE(18:0)', 'LPE(22:5)', 'PE(18:0_22:6)', 'PE(18:1_20:5)', 'PE(18:1_22:6)', 'PE(P-18:1_22:6)',
                     'LPE(18:1)', 'PE(16:0_20:5)', 'PE(16:0_22:6)', 'PE(16:1_22:6)', 'PE(18:0_20:3)', 'PE(18:0_20:4)', 'PE(18:0_20:5)')
selected_lipids <- c('PI(14:1_18:1)', 'PI(16:0_16:1)', 'PI(16:1_18:1)', 'PI(18:0_18:1)', 'PI(18:0_18:2)', 'PI(18:1_18:1)',
                     'PI(18:1_18:2)', 'PI(20:4_20:4)', 'PI(O-18:0_18:1)', 'PI(16:0_22:6)', 'PI(18:1_22:6)', 'PI(18:0_22:6)')
selected_lipids <- c('PS(16:0_16:1)', 'PS(16:0_18:1)', 'PS(16:1_18:0)', 'PS(16:1_18:1)', 'PS(18:0_18:4)', 'PS(18:1_18:1)', 'PS(18:1_18:2)', 'PS(18:1_20:4)')


#Glycerolipids
selected_lipids <- c('DG(14:0_16:0)', 'DG(14:0_18:1)', 'DG(16:0_18:0)', 'DG(16:1_18:1)', 'DG(18:0_18:0)', 'DG(16:0_18:0)', 'DG(16:1_18:1)')
selected_lipids <- c('TG(12:0_12:0_14:0)', 'TG(12:0_12:0_16:0)', 'TG(13:0_15:0_15:0)', 'TG(14:0_14:0_15:1)', 'TG(14:0_15:0_16:0)', 'TG(14:0_15:1_16:1)', 
                     'TG(15:0_16:0_16:0)', 'TG(16:0_20:3_22:6)', 'TG(16:1_16:1_20:4)', 'TG(17:0_18:0_18:0)', 'TG(18:1_20:3_22:6)',
                     'TG(16:0_18:0_22:1)', 'TG(16:0_18:0_18:0)', 'TG(16:1_18:1_18:2)', 'TG(18:0_18:1_19:1)', 'TG(16:0_16:0_18:0)', 'TG(16:0_17:0_18:1)', 'TG(16:0_18:0_20:1)',
                     'TG(18:0_18:1_20:3)', 'TG(14:0_16:0_16:0)', 'TG(14:0_16:0_18:1)', 'TG(15:0_16:0_18:1)', 'TG(16:0_16:0_16:0)',
                     'TG(16:0_16:1_20:4)', 'TG(16:0_17:1_18:1)', 'TG(16:0_18:0_22:6)', 'TG(16:0_18:1_18:1)', 'TG(16:0_18:1_18:3)')

#Sterols
selected_lipids <- c('CE(20:1)', 'CE(22:0)', 'CE(22:1)', 'CE(22:2)', 'CE(22:3)', 'CE(22:4)', 'CE(22:5)', 'CE(24:1)', 'CE(26:0)',
                     'CE(20:0)', 'CE(26:0)', 
                     'CE(14:0)', 'CE(20:3)',
                     'Cholesterol', 'Cholesterol(d7)')
#Sphingolipids
selected_lipids <- c('Hex2Cer(18:0;O2/16:0)', 'Hex2Cer(18:0;O2/22:0)', 'HexCer(18:0;O2/16:0)', 'HexCer(18:0;O2/22:0)', 'HexCer(18:0;O2/24:0)', 'HexCer(18:0;O2/24:1)',
                     'HexCer(18:1;O2/14:0)', 'HexCer(18:1;O2/24:0)', 'HexCer(36:0;O2)', 'HexCer(38:0;O2)', 'HexCer(40:1;O2)', 'HexCer(44:1;O2)', 'HexHexNAcCer(18:1;O2/24:0)',
                     'Hex2Cer(18:1;O2/18:0)', 'Hex2Cer(18:1;O2/20:0)', 'Hex2Cer(18:1;O2/24:0)', 'HexCer(18:0;O2/16:0)', 'HexCer(18:0;O2/22:0)', 'HexCer(18:0;O2/24:0)', 'HexCer(18:0;O2/24:1)',
                     'HexCer(18:1;O2/16:0)', 'HexCer(18:1;O2/22:0)', 'HexCer(18:1;O2/24:0)', 'HexCer(18:1;O2/24:1)', 'HexCer(18:1;O2/26:1)', 'HexCer(36:1;O2)', 'HexHexNAcCer(18:1;O2/16:0)', 'HexHexNAcCer(18:1;O2/22:0)', 'HexHexNAcCer(18:1;O2/24:1)')
selected_lipids <- c('SM(34:2;O2)', 'SM(36:1;O2)', 'SM(33:1;O2)', 'SM(34:1;O2)')
selected_lipids <- c('Cer(18:0;O/24:0)', 'Cer(18:0;O2/16:0)', 'Cer(18:0;O2/18:0)', 'Cer(18:0;O2/20:0)', 'Cer(18:0;O2/22:0)', 'Cer(18:0;O2/24:0)', 'Cer(18:1;O2/24:0)', 'Cer(18:1;O2/26:0)', 'Cer(18:1;O2/26:1)', 'Cer(18:2;O2/24:0)')
selected_lipids <- c('Cer(18:0;O2/16:0)', 'Cer(18:0;O2/18:0)', 'Cer(18:0;O2/20:0)', 'Cer(18:0;O2/22:0)', 'Cer(18:0;O2/24:0)', 'Cer(18:1;O2/24:0)', 'Cer(18:1;O2/26:0)', 'Cer(18:1;O2/26:1)')

#Fatty acids
selected_lipids <- c('FA(18:0)', 'FA(18:1)', 'FA(22:5)', 'FA(20:1)', 'FA(26:1)', 'CAR(18:0)', 'CAR(18:1)')


#Update to filtered data (skip if not filtering)
filtered_data <- data_matrix[ ,colnames(data_matrix) %in% selected_lipids, ]
# Update the norm component in the mSet object
mSet$dataSet$norm <- filtered_data

#Plot heatmap
mSet<-PlotStaticHeatMap(mSet, "heatmap_subset", "png", 300, width=NA, "norm", "row", "euclidean", "ward.D","bwm", 8,8, "overview", F, T, NULL, T, F, T, T, T)
