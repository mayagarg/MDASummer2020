library(Seurat)
library("devtools")
library(Matrix)
library(dplyr)

#load in the data
MUV41.data <- read.table("/Users/mayagarg/Downloads/GSE119926_RAW copy/GSM3905417_MUV41.txt", header=TRUE)

MUV41mg.data <- matrix() 
MUV41mg.data <- rbind(MUV41.data["SPARC", ], MUV41.data["P2RY12", ], MUV41.data["TMEM119", ], MUV41.data["SLC2A5", ], MUV41.data["SALL1", ], MUV41.data["ECSCR", ], 
                      MUV41.data["SCAMP5", ], MUV41.data["PTPRM", ], MUV41.data["KCND1", ], MUV41.data["SDK1", ], MUV41.data["CLSTN1", ], MUV41.data["SGCE", ], MUV41.data["SLC12A2", ],
                      MUV41.data["CADM1", ], MUV41.data["SEMA4G", ], MUV41.data["PROS1", ], MUV41.data["SLCO2B1", ], MUV41.data["RTN4RL1", ], MUV41.data["CMTM4", ], 
                      MUV41.data["FAT3", ], MUV41.data["SMO", ], MUV41.data["MRC2", ], MUV41.data["JAM2", ], MUV41.data["PLXNA4", ], MUV41.data["SLC46A1", ], 
                      MUV41.data["AGMO", ], MUV41.data["TMEM204", ], MUV41.data["BST2", ], MUV41.data["C2", ], MUV41.data["C4B", ], MUV41.data["CD74", ], 
                      MUV41.data["CFB", ], MUV41.data["CIITA", ], MUV41.data["CTSC", ], MUV41.data["EXO1", ], MUV41.data["GBP2", ], 
                      MUV41.data["IFNB1", ], MUV41.data["IRF1", ], MUV41.data["ITK", ], MUV41.data["JAK3", ], MUV41.data["LY86", ], MUV41.data["MASP1", ], 
                      MUV41.data["OAS2", ], MUV41.data["P2RY14", ], MUV41.data["SEMA4A", ], MUV41.data["SERPING1", ], MUV41.data["STAT1", ], 
                      MUV41.data["TLR2", ], MUV41.data["TLR3", ], MUV41.data["TNF", ], MUV41.data["TNFSF10", ], MUV41.data["TNFSF13b", ], MUV41.data["TNFSF15", ],MUV41.data["TNFSF8", ], 
                      MUV41.data["CX3CR1", ])


MUV41mg.data[is.na(MUV41mg.data)] <- 0

#initialize the Seurat data with the raw (non-normalized) data
#min.cells= include features (genes) detected in at least X many cells 
#min.features = include cells where at least Y many features (genes) are detected

MUV41 <- CreateSeuratObject(counts = MUV41mg.data, project= "GSM3905417_MUV41")
MUV41
#An object of class Seurat 
#55 features (genes) across 338 samples (cells) within 1 assay 
#Active assay: RNA (55 features, 0 variable features)

#QC pre processing
#MT is mitochondrial gene marker 
MUV41[["percent.mt"]] <- PercentageFeatureSet(MUV41, pattern = "^MT-")

#Normalizing data 
MUV41 <- NormalizeData(MUV41, normalization.method = "LogNormalize", scale.factor = 10000)


#identification of highly variable features (feature selection)
MUV41 <- FindVariableFeatures(MUV41, selection.method = "vst")


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(MUV41), 10)

# plot variable features with and without labels
#check to see if there are any immune related genes here
plot1 <- VariableFeaturePlot(MUV41)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
CombinePlots(plots = list(plot1, plot2))

#scaling the data
all.genes <- rownames(MUV41)
MUV41 <- ScaleData(MUV41, features = all.genes)

#perform PCA on scaled data using linear dimension reduction
#heat maps require PCA
MUV41 <- RunPCA(MUV41, features = VariableFeatures(object = MUV41))

#clustering cells
MUV41 <- FindNeighbors(MUV41, dims = 1:15)
MUV41 <- FindClusters(MUV41)

#determine marker groups and cluster data 
MUV41.markers <- FindAllMarkers(MUV41, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MUV41.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#find top markers and plot
top10 <- MUV41.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(MUV41) 


#want to extract the different clusters as individual matrices 
table(Idents(MUV41)) #how many cells are in each cluster 

# 0   1   2   3  4   5   --> cluster
# 83  66  61 59 45  24  --> number of cells per cluster

#extract cluster 0
cluster0.matrix.MUV41 <- as.matrix(GetAssayData(MUV41, slot = "counts")[, WhichCells(MUV41, ident = "0")])
ncol(cluster0.matrix.MUV41) #83 --> correct (column is cell)

#extract cluster 1
cluster1.matrix.MUV41 <- as.matrix(GetAssayData(MUV41, slot = "counts")[, WhichCells(MUV41, ident = "1")])
ncol(cluster1.matrix.MUV41) #66 --> correct  (column is cell)

#extract cluster 2
cluster2.matrix.MUV41 <- as.matrix(GetAssayData(MUV41, slot = "counts")[, WhichCells(MUV41, ident = "2")])
ncol(cluster2.matrix.MUV41) #61 --> correct  (column is cell)

#extract cluster 3
cluster3.matrix.MUV41 <- as.matrix(GetAssayData(MUV41, slot = "counts")[, WhichCells(MUV41, ident = "3")])
ncol(cluster3.matrix.MUV41) #59 --> correct  (column is cell)

#extract cluster 4
cluster4.matrix.MUV41 <- as.matrix(GetAssayData(MUV41, slot = "counts")[, WhichCells(MUV41, ident = "4")])
ncol(cluster4.matrix.MUV41) #45 --> correct  (column is cell)

#extract cluster 5
cluster5.matrix.MUV41 <- as.matrix(GetAssayData(MUV41, slot = "counts")[, WhichCells(MUV41, ident = "5")])
ncol(cluster5.matrix.MUV41) #24 --> correct  (column is cell)

#find the average expression value per cluster (sums all of the expression levels in the 
#cluster and divides by the number of cells --> average sum of expression per cell)

avg.c0.MUV41 <- sum(cluster0.matrix.MUV41)/ncol(cluster0.matrix.MUV41)
avg.c1.MUV41 <- sum(cluster1.matrix.MUV41)/ncol(cluster1.matrix.MUV41)
avg.c2.MUV41 <- sum(cluster2.matrix.MUV41)/ncol(cluster2.matrix.MUV41)
avg.c3.MUV41 <- sum(cluster3.matrix.MUV41)/ncol(cluster3.matrix.MUV41)
avg.c4.MUV41 <- sum(cluster4.matrix.MUV41)/ncol(cluster4.matrix.MUV41)
avg.c5.MUV41 <- sum(cluster5.matrix.MUV41)/ncol(cluster5.matrix.MUV41)


c("C0", "C1", "C2", "C3", "C4", "C5")[which.max(c(avg.c0.MUV41, avg.c1.MUV41, avg.c2.MUV41, avg.c3.MUV41, avg.c4.MUV41, avg.c5.MUV41))] #C3

