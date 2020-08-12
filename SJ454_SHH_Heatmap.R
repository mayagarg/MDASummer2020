library(Seurat)
library("devtools")
library(Matrix)
library(dplyr)

#load in the data
SJ454.data <- read.table("/Users/mayagarg/Downloads/GSE119926_RAW copy/GSM3905423_SJ454.txt", header=TRUE)

SJ454mg.data <- matrix() 
SJ454mg.data <- rbind(SJ454.data["SPARC", ], SJ454.data["P2RY12", ], SJ454.data["TMEM119", ], SJ454.data["SLC2A5", ], SJ454.data["SALL1", ], SJ454.data["ECSCR", ], 
                      SJ454.data["SCAMP5", ], SJ454.data["PTPRM", ], SJ454.data["KCND1", ], SJ454.data["SDK1", ], SJ454.data["CLSTN1", ], SJ454.data["SGCE", ], SJ454.data["SLC12A2", ],
                      SJ454.data["CADM1", ], SJ454.data["SEMA4G", ], SJ454.data["PROS1", ], SJ454.data["SLCO2B1", ], SJ454.data["RTN4RL1", ], SJ454.data["CMTM4", ], 
                      SJ454.data["FAT3", ], SJ454.data["SMO", ], SJ454.data["MRC2", ], SJ454.data["JAM2", ], SJ454.data["PLXNA4", ], SJ454.data["SLC46A1", ], 
                      SJ454.data["AGMO", ], SJ454.data["TMEM204", ], SJ454.data["BST2", ], SJ454.data["C2", ], SJ454.data["C4B", ], SJ454.data["CD74", ], 
                      SJ454.data["CFB", ], SJ454.data["CIITA", ], SJ454.data["CTSC", ], SJ454.data["EXO1", ], SJ454.data["GBP2", ], 
                      SJ454.data["IFNB1", ], SJ454.data["IRF1", ], SJ454.data["ITK", ], SJ454.data["JAK3", ], SJ454.data["LY86", ], SJ454.data["MASP1", ], 
                      SJ454.data["OAS2", ], SJ454.data["P2RY14", ], SJ454.data["SEMA4A", ], SJ454.data["SERPING1", ], SJ454.data["STAT1", ], 
                      SJ454.data["TLR2", ], SJ454.data["TLR3", ], SJ454.data["TNF", ], SJ454.data["TNFSF10", ], SJ454.data["TNFSF13b", ], SJ454.data["TNFSF15", ],SJ454.data["TNFSF8", ], 
                      SJ454.data["CX3CR1", ])

SJ454mg.data[is.na(SJ454mg.data)] <- 0

#initialize the Seurat data with the raw (non-normalized) data
#min.cells= include features (genes) detected in at least X many cells 
#min.features = include cells where at least Y many features (genes) are detected

SJ454 <- CreateSeuratObject(counts = SJ454mg.data, project= "GSM3905423_SJ454")
SJ454
#An object of class Seurat 
#55 features (genes) across 293 samples (cells) within 1 assay 
#Active assay: RNA (55 features, 47 variable features)

#QC pre processing
#MT is mitochondrial gene marker 
SJ454[["percent.mt"]] <- PercentageFeatureSet(SJ454, pattern = "^MT-")

#Normalizing data 
SJ454 <- NormalizeData(SJ454, normalization.method = "LogNormalize", scale.factor = 10000)

#identification of highly variable features (feature selection)
SJ454 <- FindVariableFeatures(SJ454, selection.method = "vst")

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(SJ454), 10)

# plot variable features with and without labels
#check to see if there are any immune related genes here
plot1 <- VariableFeaturePlot(SJ454)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
CombinePlots(plots = list(plot1, plot2))

#scaling the data
all.genes <- rownames(SJ454)
SJ454 <- ScaleData(SJ454, features = all.genes)

#perform PCA on scaled data using linear dimension reduction
#heat maps require PCA
SJ454 <- RunPCA(SJ454, features = VariableFeatures(object = SJ454), approx=FALSE)

#clustering cells
SJ454 <- FindNeighbors(SJ454, dims = 1:15)
SJ454 <- FindClusters(SJ454)


#determine marker groups and cluster data 
SJ454.markers <- FindAllMarkers(SJ454, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SJ454.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#find top markers and plot
top10 <- SJ454.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(SJ454) 

#want to extract the different clusters as individual matrices 
table(Idents(SJ454)) #how many cells are in each cluster 

# 0   1   2   3      --> cluster
# 101 76  64  52    --> number of cells per cluster

#extract cluster 0
cluster0.matrix.SJ454 <- as.matrix(GetAssayData(SJ454, slot = "counts")[, WhichCells(SJ454, ident = "0")])
ncol(cluster0.matrix.SJ454) #98 --> correct (column is cell)

#extract cluster 1
cluster1.matrix.SJ454 <- as.matrix(GetAssayData(SJ454, slot = "counts")[, WhichCells(SJ454, ident = "1")])
ncol(cluster1.matrix.SJ454) #81 --> correct  (column is cell)

#extract cluster 2
cluster2.matrix.SJ454 <- as.matrix(GetAssayData(SJ454, slot = "counts")[, WhichCells(SJ454, ident = "2")])
ncol(cluster2.matrix.SJ454) #66 --> correct  (column is cell)

#extract cluster 3
cluster3.matrix.SJ454 <- as.matrix(GetAssayData(SJ454, slot = "counts")[, WhichCells(SJ454, ident = "3")])
ncol(cluster3.matrix.SJ454) #48 --> correct  (column is cell)

#find the average expression value per cluster (sums all of the expression levels in the 
#cluster and divides by the number of cells --> average sum of expression per cell)

avg.c0.SJ454 <- sum(cluster0.matrix.SJ454)/ncol(cluster0.matrix.SJ454)
avg.c1.SJ454 <- sum(cluster1.matrix.SJ454)/ncol(cluster1.matrix.SJ454)
avg.c2.SJ454 <- sum(cluster2.matrix.SJ454)/ncol(cluster2.matrix.SJ454)
avg.c3.SJ454 <- sum(cluster3.matrix.SJ454)/ncol(cluster3.matrix.SJ454)


c("C0", "C1", "C2", "C3")[which.max(c(avg.c0.SJ454, avg.c1.SJ454, avg.c2.SJ454, avg.c3.SJ454))] #C2


