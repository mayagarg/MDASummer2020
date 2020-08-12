library(Seurat)
library("devtools")
library(Matrix)
library(dplyr)

#load in the data
BCH807.data <- read.table("/Users/mayagarg/Downloads/GSE119926_RAW copy/GSM3905406_BCH807.txt", header=TRUE)

BCH807mg.data <- matrix() 
BCH807mg.data <- rbind(BCH807.data["SPARC", ], BCH807.data["P2RY12", ], BCH807.data["TMEM119", ], BCH807.data["SLC2A5", ], BCH807.data["SALL1", ], BCH807.data["ECSCR", ], 
                      BCH807.data["SCAMP5", ], BCH807.data["PTPRM", ], BCH807.data["KCND1", ], BCH807.data["SDK1", ], BCH807.data["CLSTN1", ], BCH807.data["SGCE", ], BCH807.data["SLC12A2", ],
                      BCH807.data["CADM1", ], BCH807.data["SEMA4G", ], BCH807.data["PROS1", ], BCH807.data["SLCO2B1", ], BCH807.data["RTN4RL1", ], BCH807.data["CMTM4", ], 
                      BCH807.data["FAT3", ], BCH807.data["SMO", ], BCH807.data["MRC2", ], BCH807.data["JAM2", ], BCH807.data["PLXNA4", ], BCH807.data["SLC46A1", ], 
                      BCH807.data["AGMO", ], BCH807.data["TMEM204", ], BCH807.data["BST2", ], BCH807.data["C2", ], BCH807.data["C4B", ], BCH807.data["CD74", ], 
                      BCH807.data["CFB", ], BCH807.data["CIITA", ], BCH807.data["CTSC", ], BCH807.data["EXO1", ], BCH807.data["GBP2", ], 
                      BCH807.data["IFNB1", ], BCH807.data["IRF1", ], BCH807.data["ITK", ], BCH807.data["JAK3", ], BCH807.data["LY86", ], BCH807.data["MASP1", ], 
                      BCH807.data["OAS2", ], BCH807.data["P2RY14", ], BCH807.data["SEMA4A", ], BCH807.data["SERPING1", ], BCH807.data["STAT1", ], 
                      BCH807.data["TLR2", ], BCH807.data["TLR3", ], BCH807.data["TNF", ], BCH807.data["TNFSF10", ], BCH807.data["TNFSF13b", ], BCH807.data["TNFSF15", ],BCH807.data["TNFSF8", ], 
                      BCH807.data["CX3CR1", ]) 


BCH807mg.data[is.na(BCH807mg.data)] <- 0

#initialize the Seurat data with the raw (non-normalized) data
#min.cells= include features (genes) detected in at least X many cells 
#min.features = include cells where at least Y many features (genes) are detected

BCH807 <- CreateSeuratObject(counts = BCH807mg.data, project= "GSM3905406_BCH807")
BCH807
#An object of class Seurat 
#55 features (genes) across 307 samples (cells) within 1 assay 
#Active assay: RNA (55 features, 0 variable features)

#QC pre processing
#MT is mitochondrial gene marker 
BCH807[["percent.mt"]] <- PercentageFeatureSet(BCH807, pattern = "^MT-")

#Normalizing data 
BCH807 <- NormalizeData(BCH807, normalization.method = "LogNormalize", scale.factor = 10000)


#identification of highly variable features (feature selection)
BCH807 <- FindVariableFeatures(BCH807, selection.method = "vst")


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(BCH807), 10)

# plot variable features with and without labels
#check to see if there are any immune related genes here
plot1 <- VariableFeaturePlot(BCH807)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
CombinePlots(plots = list(plot1, plot2))

#scaling the data
all.genes <- rownames(BCH807)
BCH807 <- ScaleData(BCH807, features = all.genes)

#perform PCA on scaled data using linear dimension reduction
#heat maps require PCA
BCH807 <- RunPCA(BCH807, features = VariableFeatures(object = BCH807), approx=FALSE)


#clustering cells
BCH807 <- FindNeighbors(BCH807, dims = 1:15)
BCH807 <- FindClusters(BCH807)

#determine marker groups and cluster data 
BCH807.markers <- FindAllMarkers(BCH807, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
BCH807.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)



#find top markers and plot
top10 <- BCH807.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(BCH807) #+ theme(axis.text.y = element_text(size=3))

#want to extract the different clusters as individual matrices 
table(Idents(BCH807)) #how many cells are in each cluster 

# 0  1  2  3 --> cluster
# 101 87 77 42 --> number of cells per cluster

#extract cluster 0
cluster0.matrix.BCH807 <- as.matrix(GetAssayData(BCH807, slot = "counts")[, WhichCells(BCH807, ident = "0")])
ncol(cluster0.matrix.BCH807) #101 --> correct (column is cell)

#extract cluster 1
cluster1.matrix.BCH807 <- as.matrix(GetAssayData(BCH807, slot = "counts")[, WhichCells(BCH807, ident = "1")])
ncol(cluster1.matrix.BCH807) #87 --> correct  (column is cell)

#extract cluster 2
cluster2.matrix.BCH807 <- as.matrix(GetAssayData(BCH807, slot = "counts")[, WhichCells(BCH807, ident = "2")])
ncol(cluster2.matrix.BCH807) #77 --> correct  (column is cell)

#extract cluster 3 
cluster3.matrix.BCH807 <- as.matrix(GetAssayData(BCH807, slot = "counts")[, WhichCells(BCH807, ident = "3")])
ncol(cluster3.matrix.BCH807) #42 --> correct  (column is cell)

#find the avergae expression value per cluster (sums all of the expression levels in the 
#cluster and divides by the number of cells --> average sum of expression per cell)

avg.c0.BCH807 <- sum(cluster0.matrix.BCH807)/ncol(cluster0.matrix.BCH807)
avg.c1.BCH807 <- sum(cluster1.matrix.BCH807)/ncol(cluster1.matrix.BCH807)
avg.c2.BCH807 <- sum(cluster2.matrix.BCH807)/ncol(cluster2.matrix.BCH807)
avg.c3.BCH807 <- sum(cluster3.matrix.BCH807)/ncol(cluster3.matrix.BCH807)

c("C0", "C1", "C2", "C3")[which.max(c(avg.c0.BCH807, avg.c1.BCH807, avg.c2.BCH807, avg.c3.BCH807))] #C1

