library(Seurat)
library("devtools")
library(Matrix)
library(dplyr)

#load in the data
MUV29.data <- read.table("/Users/mayagarg/Downloads/GSE119926_RAW copy/GSM3905413_MUV29.txt", header=TRUE)

MUV29mg.data <- matrix() 
MUV29mg.data <- rbind(MUV29.data["SPARC", ], MUV29.data["P2RY12", ], MUV29.data["TMEM119", ], MUV29.data["SLC2A5", ], MUV29.data["SALL1", ], MUV29.data["ECSCR", ], 
                      MUV29.data["SCAMP5", ], MUV29.data["PTPRM", ], MUV29.data["KCND1", ], MUV29.data["SDK1", ], MUV29.data["CLSTN1", ], MUV29.data["SGCE", ], MUV29.data["SLC12A2", ],
                      MUV29.data["CADM1", ], MUV29.data["SEMA4G", ], MUV29.data["PROS1", ], MUV29.data["SLCO2B1", ], MUV29.data["RTN4RL1", ], MUV29.data["CMTM4", ], 
                      MUV29.data["FAT3", ], MUV29.data["SMO", ], MUV29.data["MRC2", ], MUV29.data["JAM2", ], MUV29.data["PLXNA4", ], MUV29.data["SLC46A1", ], 
                      MUV29.data["AGMO", ], MUV29.data["TMEM204", ], MUV29.data["BST2", ], MUV29.data["C2", ], MUV29.data["C4B", ], MUV29.data["CD74", ], 
                      MUV29.data["CFB", ], MUV29.data["CIITA", ], MUV29.data["CTSC", ], MUV29.data["EXO1", ], MUV29.data["GBP2", ], 
                      MUV29.data["IFNB1", ], MUV29.data["IRF1", ], MUV29.data["ITK", ], MUV29.data["JAK3", ], MUV29.data["LY86", ], MUV29.data["MASP1", ], 
                      MUV29.data["OAS2", ], MUV29.data["P2RY14", ], MUV29.data["SEMA4A", ], MUV29.data["SERPING1", ], MUV29.data["STAT1", ], 
                      MUV29.data["TLR2", ], MUV29.data["TLR3", ], MUV29.data["TNF", ], MUV29.data["TNFSF10", ], MUV29.data["TNFSF13b", ], MUV29.data["TNFSF15", ],MUV29.data["TNFSF8", ], 
                      MUV29.data["CX3CR1", ]) 


MUV29mg.data[is.na(MUV29mg.data)] <- 0

#initialize the Seurat data with the raw (non-normalized) data
#min.cells= include features (genes) detected in at least X many cells 
#min.features = include cells where at least Y many features (genes) are detected

MUV29 <- CreateSeuratObject(counts = MUV29mg.data, project= "GSM3905413_MUV29")
MUV29
#An object of class Seurat 
#55 features (genes) across 314 samples (cells) within 1 assay 
#Active assay: RNA (55 features, 0 variable features)

#QC pre processing
#MT is mitochondrial gene marker 
MUV29[["percent.mt"]] <- PercentageFeatureSet(MUV29, pattern = "^MT-")

#Normalizing data 
MUV29 <- NormalizeData(MUV29, normalization.method = "LogNormalize", scale.factor = 10000)


#identification of highly variable features (feature selection)
MUV29 <- FindVariableFeatures(MUV29, selection.method = "vst")


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(MUV29), 10)

# plot variable features with and without labels
#check to see if there are any immune related genes here
plot1 <- VariableFeaturePlot(MUV29)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
CombinePlots(plots = list(plot1, plot2))

#scaling the data
all.genes <- rownames(MUV29)
MUV29 <- ScaleData(MUV29, features = all.genes)

#perform PCA on scaled data using linear dimension reduction
#heat maps require PCA
MUV29 <- RunPCA(MUV29, features = VariableFeatures(object = MUV29), approx=FALSE)


#clustering cells
MUV29 <- FindNeighbors(MUV29, dims = 1:15)
MUV29 <- FindClusters(MUV29)

#determine marker groups and cluster data 
MUV29.markers <- FindAllMarkers(MUV29, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MUV29.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)


#find top markers and plot
top10 <- MUV29.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(MUV29) #+ theme(axis.text.y = element_text(size=3))

#want to extract the different clusters as individual matrices 
table(Idents(MUV29)) #how many cells are in each cluster 

# 0   1   2   --> cluster
# 126 112 76 --> number of cells per cluster

#extract cluster 0
cluster0.matrix.MUV29 <- as.matrix(GetAssayData(MUV29, slot = "counts")[, WhichCells(MUV29, ident = "0")])
ncol(cluster0.matrix.MUV29) #126 --> correct (column is cell)

#extract cluster 1
cluster1.matrix.MUV29 <- as.matrix(GetAssayData(MUV29, slot = "counts")[, WhichCells(MUV29, ident = "1")])
ncol(cluster1.matrix.MUV29) #112 --> correct  (column is cell)

#extract cluster 2
cluster2.matrix.MUV29 <- as.matrix(GetAssayData(MUV29, slot = "counts")[, WhichCells(MUV29, ident = "2")])
ncol(cluster2.matrix.MUV29) #76 --> correct  (column is cell)

#find the average expression value per cluster (sums all of the expression levels in the 
#cluster and divides by the number of cells --> average sum of expression per cell)

avg.c0.MUV29 <- sum(cluster0.matrix.MUV29)/ncol(cluster0.matrix.MUV29)
avg.c1.MUV29 <- sum(cluster1.matrix.MUV29)/ncol(cluster1.matrix.MUV29)
avg.c2.MUV29 <- sum(cluster2.matrix.MUV29)/ncol(cluster2.matrix.MUV29)

c("C0", "C1", "C2", "C3")[which.max(c(avg.c0.MUV29, avg.c1.MUV29, avg.c2.MUV29))] #C0


