library(Seurat)
library("devtools")
library(Matrix)
library(dplyr)

#load in the data
BCH1031.data <- read.table("/Users/mayagarg/Downloads/GSE119926_RAW copy/GSM3905408_BCH1031.txt", header=TRUE)

BCH1031mg.data <- matrix() 
BCH1031mg.data <- rbind(BCH1031.data["SPARC", ], BCH1031.data["P2RY12", ], BCH1031.data["TMEM119", ], BCH1031.data["SLC2A5", ], BCH1031.data["SALL1", ], BCH1031.data["ECSCR", ], 
                        BCH1031.data["SCAMP5", ], BCH1031.data["PTPRM", ], BCH1031.data["KCND1", ], BCH1031.data["SDK1", ], BCH1031.data["CLSTN1", ], BCH1031.data["SGCE", ], BCH1031.data["SLC12A2", ],
                        BCH1031.data["CADM1", ], BCH1031.data["SEMA4G", ], BCH1031.data["PROS1", ], BCH1031.data["SLCO2B1", ], BCH1031.data["RTN4RL1", ], BCH1031.data["CMTM4", ], 
                        BCH1031.data["FAT3", ], BCH1031.data["SMO", ], BCH1031.data["MRC2", ], BCH1031.data["JAM2", ], BCH1031.data["PLXNA4", ], BCH1031.data["SLC46A1", ], 
                        BCH1031.data["AGMO", ], BCH1031.data["TMEM204", ], BCH1031.data["BST2", ], BCH1031.data["C2", ], BCH1031.data["C4B", ], BCH1031.data["CD74", ], 
                        BCH1031.data["CFB", ], BCH1031.data["CIITA", ], BCH1031.data["CTSC", ], BCH1031.data["EXO1", ], BCH1031.data["GBP2", ], 
                        BCH1031.data["IFNB1", ], BCH1031.data["IRF1", ], BCH1031.data["ITK", ], BCH1031.data["JAK3", ], BCH1031.data["LY86", ], BCH1031.data["MASP1", ], 
                        BCH1031.data["OAS2", ], BCH1031.data["P2RY14", ], BCH1031.data["SEMA4A", ], BCH1031.data["SERPING1", ], BCH1031.data["STAT1", ], 
                        BCH1031.data["TLR2", ], BCH1031.data["TLR3", ], BCH1031.data["TNF", ], BCH1031.data["TNFSF10", ], BCH1031.data["TNFSF13b", ], BCH1031.data["TNFSF15", ],BCH1031.data["TNFSF8", ], 
                        BCH1031.data["CX3CR1", ]) 
BCH1031mg.data[is.na(BCH1031mg.data)] <- 0

#initialize the Seurat data with the raw (non-normalized) data
#min.cells= include features (genes) detected in at least X many cells 
#min.features = include cells where at least Y many features (genes) are detected

BCH1031 <- CreateSeuratObject(counts = BCH1031mg.data, project= "GSM3905408_BCH1031")
BCH1031
#An object of class Seurat 
#55 features (genes) across 273 samples (cells) within 1 assay 
#Active assay: RNA (54 features, 0 variable features)

#QC pre processing
#MT is mitochondrial gene marker 
BCH1031[["percent.mt"]] <- PercentageFeatureSet(BCH1031, pattern = "^MT-")

#Normalizing data 
BCH1031 <- NormalizeData(BCH1031, normalization.method = "LogNormalize", scale.factor = 10000)


#identification of highly variable features (feature selection)
BCH1031 <- FindVariableFeatures(BCH1031, selection.method = "vst")


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(BCH1031), 10)

# plot variable features with and without labels
#check to see if there are any immune related genes here
plot1 <- VariableFeaturePlot(BCH1031)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
CombinePlots(plots = list(plot1, plot2))

#scaling the data
all.genes <- rownames(BCH1031)
BCH1031 <- ScaleData(BCH1031, features = all.genes)

#perform PCA on scaled data using linear dimension reduction
#heat maps require PCA
BCH1031 <- RunPCA(BCH1031, features = VariableFeatures(object = BCH1031), approx=FALSE)


#clustering cells
BCH1031 <- FindNeighbors(BCH1031, dims = 1:15)
BCH1031 <- FindClusters(BCH1031)

#determine marker groups and cluster data 
BCH1031.markers <- FindAllMarkers(BCH1031, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
BCH1031.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#find top markers and plot
top10 <- BCH1031.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(BCH1031) #+ theme(axis.text.y = element_text(size=3))

#want to extract the different clusters as individual matrices 
table(Idents(BCH1031)) #how many cells are in each cluster 

# 0   1   2  --> cluster
# 112 86  75 --> number of cells per cluster

#extract cluster 0
cluster0.matrix.BCH1031 <- as.matrix(GetAssayData(BCH1031, slot = "counts")[, WhichCells(BCH1031, ident = "0")])
ncol(cluster0.matrix.BCH1031) #112 --> correct (column is cell)

#extract cluster 1
cluster1.matrix.BCH1031 <- as.matrix(GetAssayData(BCH1031, slot = "counts")[, WhichCells(BCH1031, ident = "1")])
ncol(cluster1.matrix.BCH1031) #86 --> correct  (column is cell)

#extract cluster 2
cluster2.matrix.BCH1031 <- as.matrix(GetAssayData(BCH1031, slot = "counts")[, WhichCells(BCH1031, ident = "2")])
ncol(cluster2.matrix.BCH1031) #95 --> correct  (column is cell)

#find the average expression value per cluster (sums all of the expression levels in the 
#cluster and divides by the number of cells --> average sum of expression per cell)

avg.c0.BCH1031 <- sum(cluster0.matrix.BCH1031)/ncol(cluster0.matrix.BCH1031)
avg.c1.BCH1031 <- sum(cluster1.matrix.BCH1031)/ncol(cluster1.matrix.BCH1031)
avg.c2.BCH1031 <- sum(cluster2.matrix.BCH1031)/ncol(cluster2.matrix.BCH1031)

c("C0", "C1", "C2", "C3")[which.max(c(avg.c0.BCH1031, avg.c1.BCH1031, avg.c2.BCH1031))] #C0



