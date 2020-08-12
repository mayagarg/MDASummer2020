library(Seurat)
library("devtools")
library(Matrix)
library(dplyr)

#load in the data
MUV19.data <- read.table("/Users/mayagarg/Downloads/GSE119926_RAW copy/GSM3905411_MUV19.txt", header=TRUE)

MUV19mg.data <- matrix() 
MUV19mg.data <- rbind(MUV19.data["SPARC", ], MUV19.data["P2RY12", ], MUV19.data["TMEM119", ], MUV19.data["SLC2A5", ], MUV19.data["SALL1", ], MUV19.data["ECSCR", ], 
                      MUV19.data["SCAMP5", ], MUV19.data["PTPRM", ], MUV19.data["KCND1", ], MUV19.data["SDK1", ], MUV19.data["CLSTN1", ], MUV19.data["SGCE", ], MUV19.data["SLC12A2", ],
                      MUV19.data["CADM1", ], MUV19.data["SEMA4G", ], MUV19.data["PROS1", ], MUV19.data["SLCO2B1", ], MUV19.data["RTN4RL1", ], MUV19.data["CMTM4", ], 
                      MUV19.data["FAT3", ], MUV19.data["SMO", ], MUV19.data["MRC2", ], MUV19.data["JAM2", ], MUV19.data["PLXNA4", ], MUV19.data["SLC46A1", ], 
                      MUV19.data["AGMO", ], MUV19.data["TMEM204", ], MUV19.data["BST2", ], MUV19.data["C2", ], MUV19.data["C4B", ], MUV19.data["CD74", ], 
                      MUV19.data["CFB", ], MUV19.data["CIITA", ], MUV19.data["CTSC", ], MUV19.data["EXO1", ], MUV19.data["GBP2", ], 
                      MUV19.data["IFNB1", ], MUV19.data["IRF1", ], MUV19.data["ITK", ], MUV19.data["JAK3", ], MUV19.data["LY86", ], MUV19.data["MASP1", ], 
                      MUV19.data["OAS2", ], MUV19.data["P2RY14", ], MUV19.data["SEMA4A", ], MUV19.data["SERPING1", ], MUV19.data["STAT1", ], 
                      MUV19.data["TLR2", ], MUV19.data["TLR3", ], MUV19.data["TNF", ], MUV19.data["TNFSF10", ], MUV19.data["TNFSF13b", ], MUV19.data["TNFSF15", ],MUV19.data["TNFSF8", ], 
                      MUV19.data["CX3CR1", ]) 


MUV19mg.data[is.na(MUV19mg.data)] <- 0

#initialize the Seurat data with the raw (non-normalized) data
#min.cells= include features (genes) detected in at least X many cells 
#min.features = include cells where at least Y many features (genes) are detected

MUV19 <- CreateSeuratObject(counts = MUV19mg.data, project= "GSM3905419_MUV19")
MUV19
#An object of class Seurat 
#55 features (genes) across 304 samples (cells) within 1 assay 
#Active assay: RNA (55 features, 0 variable features)

#QC pre processing
#MT is mitochondrial gene marker 
MUV19[["percent.mt"]] <- PercentageFeatureSet(MUV19, pattern = "^MT-")

#Normalizing data 
MUV19 <- NormalizeData(MUV19, normalization.method = "LogNormalize", scale.factor = 10000)


#identification of highly variable features (feature selection)
MUV19 <- FindVariableFeatures(MUV19, selection.method = "vst")


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(MUV19), 10)

# plot variable features with and without labels
#check to see if there are any immune related genes here
plot1 <- VariableFeaturePlot(MUV19)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
CombinePlots(plots = list(plot1, plot2))

#scaling the data
all.genes <- rownames(MUV19)
MUV19 <- ScaleData(MUV19, features = all.genes)

#perform PCA on scaled data using linear dimension reduction
#heat maps require PCA
MUV19 <- RunPCA(MUV19, features = VariableFeatures(object = MUV19), approx=FALSE)


#clustering cells
MUV19 <- FindNeighbors(MUV19, dims = 1:15)
MUV19 <- FindClusters(MUV19)

#determine marker groups and cluster data 
MUV19.markers <- FindAllMarkers(MUV19, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MUV19.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#find top markers and plot
top10 <- MUV19.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(MUV19) #+ theme(axis.text.y = element_text(size=3))

#want to extract the different clusters as individual matrices 
table(Idents(MUV19)) #how many cells are in each cluster 

# 0   1   2   --> cluster
# 126 106 72  --> number of cells per cluster

#extract cluster 0
cluster0.matrix.MUV19 <- as.matrix(GetAssayData(MUV19, slot = "counts")[, WhichCells(MUV19, ident = "0")])
ncol(cluster0.matrix.MUV19) #115 --> correct (column is cell)

#extract cluster 1
cluster1.matrix.MUV19 <- as.matrix(GetAssayData(MUV19, slot = "counts")[, WhichCells(MUV19, ident = "1")])
ncol(cluster1.matrix.MUV19) #101 --> correct  (column is cell)

#extract cluster 2
cluster2.matrix.MUV19 <- as.matrix(GetAssayData(MUV19, slot = "counts")[, WhichCells(MUV19, ident = "2")])
ncol(cluster2.matrix.MUV19) #88 --> correct  (column is cell)

#find the average expression value per cluster (sums all of the expression levels in the 
#cluster and divides by the number of cells --> average sum of expression per cell)

avg.c0.MUV19 <- sum(cluster0.matrix.MUV19)/ncol(cluster0.matrix.MUV19)
avg.c1.MUV19 <- sum(cluster1.matrix.MUV19)/ncol(cluster1.matrix.MUV19)
avg.c2.MUV19 <- sum(cluster2.matrix.MUV19)/ncol(cluster2.matrix.MUV19)

c("C0", "C1", "C2", "C3")[which.max(c(avg.c0.MUV19, avg.c1.MUV19, avg.c2.MUV19))] #C0




