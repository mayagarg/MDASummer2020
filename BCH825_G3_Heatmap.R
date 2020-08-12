library(Seurat)
library("devtools")
library(Matrix)
library(dplyr)

#load in the data
BCH825.data <- read.table("/Users/mayagarg/Downloads/GSE119926_RAW copy/GSM3905407_BCH825.txt", header=TRUE)

BCH825mg.data <- matrix() 
BCH825mg.data <- rbind(BCH825.data["SPARC", ], BCH825.data["P2RY12", ], BCH825.data["TMEM119", ], BCH825.data["SLC2A5", ], BCH825.data["SALL1", ], BCH825.data["ECSCR", ], 
                       BCH825.data["SCAMP5", ], BCH825.data["PTPRM", ], BCH825.data["KCND1", ], BCH825.data["SDK1", ], BCH825.data["CLSTN1", ], BCH825.data["SGCE", ], BCH825.data["SLC12A2", ],
                       BCH825.data["CADM1", ], BCH825.data["SEMA4G", ], BCH825.data["PROS1", ], BCH825.data["SLCO2B1", ], BCH825.data["RTN4RL1", ], BCH825.data["CMTM4", ], 
                       BCH825.data["FAT3", ], BCH825.data["SMO", ], BCH825.data["MRC2", ], BCH825.data["JAM2", ], BCH825.data["PLXNA4", ], BCH825.data["SLC46A1", ], 
                       BCH825.data["AGMO", ], BCH825.data["TMEM204", ], BCH825.data["BST2", ], BCH825.data["C2", ], BCH825.data["C4B", ], BCH825.data["CD74", ], 
                       BCH825.data["CFB", ], BCH825.data["CIITA", ], BCH825.data["CTSC", ], BCH825.data["EXO1", ], BCH825.data["GBP2", ], 
                       BCH825.data["IFNB1", ], BCH825.data["IRF1", ], BCH825.data["ITK", ], BCH825.data["JAK3", ], BCH825.data["LY86", ], BCH825.data["MASP1", ], 
                       BCH825.data["OAS2", ], BCH825.data["P2RY14", ], BCH825.data["SEMA4A", ], BCH825.data["SERPING1", ], BCH825.data["STAT1", ], 
                       BCH825.data["TLR2", ], BCH825.data["TLR3", ], BCH825.data["TNF", ], BCH825.data["TNFSF10", ], BCH825.data["TNFSF13b", ], BCH825.data["TNFSF15", ],BCH825.data["TNFSF8", ], 
                       BCH825.data["CX3CR1", ]) 

BCH825mg.data[is.na(BCH825mg.data)] <- 0

#initialize the Seurat data with the raw (non-normalized) data
#min.cells= include features (genes) detected in at least X many cells 
#min.features = include cells where at least Y many features (genes) are detected

BCH825 <- CreateSeuratObject(counts = BCH825mg.data, project= "GSM3905407_BCH825")
BCH825
#An object of class Seurat 
#55 features (genes) across 330 samples (cells) within 1 assay 
#Active assay: RNA (55 features, 0 variable features)

#QC pre processing
#MT is mitochondrial gene marker 
BCH825[["percent.mt"]] <- PercentageFeatureSet(BCH825, pattern = "^MT-")

#Normalizing data 
BCH825 <- NormalizeData(BCH825, normalization.method = "LogNormalize", scale.factor = 10000)


#identification of highly variable features (feature selection)
BCH825 <- FindVariableFeatures(BCH825, selection.method = "vst")


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(BCH825), 10)

# plot variable features with and without labels
#check to see if there are any immune related genes here
plot1 <- VariableFeaturePlot(BCH825)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
CombinePlots(plots = list(plot1, plot2))

#scaling the data
all.genes <- rownames(BCH825)
BCH825 <- ScaleData(BCH825, features = all.genes)

#perform PCA on scaled data using linear dimension reduction
#heat maps require PCA
BCH825 <- RunPCA(BCH825, features = VariableFeatures(object = BCH825), approx=FALSE)


#clustering cells
BCH825 <- FindNeighbors(BCH825, dims = 1:15)
BCH825 <- FindClusters(BCH825)

#determine marker groups and cluster data 
BCH825.markers <- FindAllMarkers(BCH825, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
BCH825.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)



#find top markers and plot
top10 <- BCH825.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(BCH825) #+ theme(axis.text.y = element_text(size=3))

#want to extract the different clusters as individual matrices 
table(Idents(BCH825)) #how many cells are in each cluster 

# 0    1   2  --> cluster
# 131 104 95 --> number of cells per cluster

#extract cluster 0
cluster0.matrix.BCH825 <- as.matrix(GetAssayData(BCH825, slot = "counts")[, WhichCells(BCH825, ident = "0")])
ncol(cluster0.matrix.BCH825) #131 --> correct (column is cell)

#extract cluster 1
cluster1.matrix.BCH825 <- as.matrix(GetAssayData(BCH825, slot = "counts")[, WhichCells(BCH825, ident = "1")])
ncol(cluster1.matrix.BCH825) #104 --> correct  (column is cell)

#extract cluster 2
cluster2.matrix.BCH825 <- as.matrix(GetAssayData(BCH825, slot = "counts")[, WhichCells(BCH825, ident = "2")])
ncol(cluster2.matrix.BCH825) #95 --> correct  (column is cell)

#find the avergae expression value per cluster (sums all of the expression levels in the 
#cluster and divides by the number of cells --> average sum of expression per cell)

avg.c0.BCH825 <- sum(cluster0.matrix.BCH825)/ncol(cluster0.matrix.BCH825)
avg.c1.BCH825 <- sum(cluster1.matrix.BCH825)/ncol(cluster1.matrix.BCH825)
avg.c2.BCH825 <- sum(cluster2.matrix.BCH825)/ncol(cluster2.matrix.BCH825)

c("C0", "C1", "C2", "C3")[which.max(c(avg.c0.BCH825, avg.c1.BCH825, avg.c2.BCH825))] #C1


