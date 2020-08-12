library(Seurat)
library("devtools")
library(Matrix)
library(dplyr)

#load in the data
SJ99.data <- read.table("/Users/mayagarg/Downloads/GSE119926_RAW copy/GSM3905420_SJ99.txt", header=TRUE)

SJ99mg.data <- matrix() 
SJ99mg.data <-  rbind(SJ99.data["SPARC", ], SJ99.data["P2RY12", ], SJ99.data["TMEM119", ], SJ99.data["SLC2A5", ], SJ99.data["SALL1", ], SJ99.data["ECSCR", ], 
                      SJ99.data["SCAMP5", ], SJ99.data["PTPRM", ], SJ99.data["KCND1", ], SJ99.data["SDK1", ], SJ99.data["CLSTN1", ], SJ99.data["SGCE", ], SJ99.data["SLC12A2", ],
                      SJ99.data["CADM1", ], SJ99.data["SEMA4G", ], SJ99.data["PROS1", ], SJ99.data["SLCO2B1", ], SJ99.data["RTN4RL1", ], SJ99.data["CMTM4", ], 
                      SJ99.data["FAT3", ], SJ99.data["SMO", ], SJ99.data["MRC2", ], SJ99.data["JAM2", ], SJ99.data["PLXNA4", ], SJ99.data["SLC46A1", ], 
                      SJ99.data["AGMO", ], SJ99.data["TMEM204", ], SJ99.data["BST2", ], SJ99.data["C2", ], SJ99.data["C4B", ], SJ99.data["CD74", ], 
                      SJ99.data["CFB", ], SJ99.data["CIITA", ], SJ99.data["CTSC", ], SJ99.data["EXO1", ], SJ99.data["GBP2", ], 
                      SJ99.data["IFNB1", ], SJ99.data["IRF1", ], SJ99.data["ITK", ], SJ99.data["JAK3", ], SJ99.data["LY86", ], SJ99.data["MASP1", ], 
                      SJ99.data["OAS2", ], SJ99.data["P2RY14", ], SJ99.data["SEMA4A", ], SJ99.data["SERPING1", ], SJ99.data["STAT1", ], 
                      SJ99.data["TLR2", ], SJ99.data["TLR3", ], SJ99.data["TNF", ], SJ99.data["TNFSF10", ], SJ99.data["TNFSF13b", ], SJ99.data["TNFSF15", ],SJ99.data["TNFSF8", ], 
                      SJ99.data["CX3CR1", ])

SJ99mg.data[is.na(SJ99mg.data)] <- 0

#initialize the Seurat data with the raw (non-normalized) data
#min.cells= include features (genes) detected in at least X many cells 
#min.features = include cells where at least Y many features (genes) are detected

SJ99 <- CreateSeuratObject(counts = SJ99mg.data, project= "GSM3905420_SJ99")
SJ99
#An object of class Seurat 
#55 features (genes) across 338 samples (cells) within 1 assay 
#Active assay: RNA (55 features, 0 variable features)

#QC pre processing
#MT is mitochondrial gene marker 
SJ99[["percent.mt"]] <- PercentageFeatureSet(SJ99, pattern = "^MT-")

#Normalizing data 
SJ99 <- NormalizeData(SJ99, normalization.method = "LogNormalize", scale.factor = 10000)

#identification of highly variable features (feature selection)
SJ99 <- FindVariableFeatures(SJ99, selection.method = "vst")

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(SJ99), 10)

# plot variable features with and without labels
#check to see if there are any immune related genes here
plot1 <- VariableFeaturePlot(SJ99)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
CombinePlots(plots = list(plot1, plot2))

#scaling the data
all.genes <- rownames(SJ99)
SJ99 <- ScaleData(SJ99, features = all.genes)

#perform PCA on scaled data using linear dimension reduction
#heat maps require PCA
SJ99 <- RunPCA(SJ99, features = VariableFeatures(object = SJ99), approx=FALSE)

#clustering cells
SJ99 <- FindNeighbors(SJ99, dims = 1:15)
SJ99 <- FindClusters(SJ99)

#determine marker groups and cluster data 
SJ99.markers <- FindAllMarkers(SJ99, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SJ99.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#find top markers and plot
top10 <- SJ99.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(SJ99) #+ theme(axis.text.y = element_text(size=3))

#want to extract the different clusters as individual matrices 
table(Idents(SJ99)) #how many cells are in each cluster 

# 0   1    2     --> cluster
# 126 110  102 --> number of cells per cluster

#extract cluster 0
cluster0.matrix.SJ99 <- as.matrix(GetAssayData(SJ99, slot = "counts")[, WhichCells(SJ99, ident = "0")])
ncol(cluster0.matrix.SJ99) #126 --> correct (column is cell)

#extract cluster 1
cluster1.matrix.SJ99 <- as.matrix(GetAssayData(SJ99, slot = "counts")[, WhichCells(SJ99, ident = "1")])
ncol(cluster1.matrix.SJ99) #110 --> correct  (column is cell)

#extract cluster 2
cluster2.matrix.SJ99 <- as.matrix(GetAssayData(SJ99, slot = "counts")[, WhichCells(SJ99, ident = "2")])
ncol(cluster2.matrix.SJ99) #102 --> correct  (column is cell)

#find the average expression value per cluster (sums all of the expression levels in the 
#cluster and divides by the number of cells --> average sum of expression per cell)

avg.c0.SJ99 <- sum(cluster0.matrix.SJ99)/ncol(cluster0.matrix.SJ99)
avg.c1.SJ99 <- sum(cluster1.matrix.SJ99)/ncol(cluster1.matrix.SJ99)
avg.c2.SJ99 <- sum(cluster2.matrix.SJ99)/ncol(cluster2.matrix.SJ99)


c("C0", "C1", "C2")[which.max(c(avg.c0.SJ99, avg.c1.SJ99, avg.c2.SJ99))] #C0



