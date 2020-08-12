library(Seurat)
library("devtools")
library(Matrix)
library(dplyr)

#load in the data
MUV27.data <- read.table("/Users/mayagarg/Downloads/GSE119926_RAW copy/GSM3905412_MUV27.txt", header=TRUE)

MUV27mg.data <- matrix() 
MUV27mg.data <- rbind(MUV27.data["SPARC", ], MUV27.data["P2RY12", ], MUV27.data["TMEM119", ], MUV27.data["SLC2A5", ], MUV27.data["SALL1", ], MUV27.data["ECSCR", ], 
                      MUV27.data["SCAMP5", ], MUV27.data["PTPRM", ], MUV27.data["KCND1", ], MUV27.data["SDK1", ], MUV27.data["CLSTN1", ], MUV27.data["SGCE", ], MUV27.data["SLC12A2", ],
                      MUV27.data["CADM1", ], MUV27.data["SEMA4G", ], MUV27.data["PROS1", ], MUV27.data["SLCO2B1", ], MUV27.data["RTN4RL1", ], MUV27.data["CMTM4", ], 
                      MUV27.data["FAT3", ], MUV27.data["SMO", ], MUV27.data["MRC2", ], MUV27.data["JAM2", ], MUV27.data["PLXNA4", ], MUV27.data["SLC46A1", ], 
                      MUV27.data["AGMO", ], MUV27.data["TMEM204", ], MUV27.data["BST2", ], MUV27.data["C2", ], MUV27.data["C4B", ], MUV27.data["CD74", ], 
                      MUV27.data["CFB", ], MUV27.data["CIITA", ], MUV27.data["CTSC", ], MUV27.data["EXO1", ], MUV27.data["GBP2", ], 
                      MUV27.data["IFNB1", ], MUV27.data["IRF1", ], MUV27.data["ITK", ], MUV27.data["JAK3", ], MUV27.data["LY86", ], MUV27.data["MASP1", ], 
                      MUV27.data["OAS2", ], MUV27.data["P2RY14", ], MUV27.data["SEMA4A", ], MUV27.data["SERPING1", ], MUV27.data["STAT1", ], 
                      MUV27.data["TLR2", ], MUV27.data["TLR3", ], MUV27.data["TNF", ], MUV27.data["TNFSF10", ], MUV27.data["TNFSF13b", ], MUV27.data["TNFSF15", ],MUV27.data["TNFSF8", ], 
                      MUV27.data["CX3CR1", ]) 


MUV27mg.data[is.na(MUV27mg.data)] <- 0

#initialize the Seurat data with the raw (non-normalized) data
#min.cells= include features (genes) detected in at least X many cells 
#min.features = include cells where at least Y many features (genes) are detected

MUV27 <- CreateSeuratObject(counts = MUV27mg.data, project= "GSM3905412_MUV27")
MUV27
#An object of class Seurat 
#55 features (genes) across 314 samples (cells) within 1 assay 
#Active assay: RNA (55 features, 0 variable features)

#QC pre processing
#MT is mitochondrial gene marker 
MUV27[["percent.mt"]] <- PercentageFeatureSet(MUV27, pattern = "^MT-")

#Normalizing data 
MUV27 <- NormalizeData(MUV27, normalization.method = "LogNormalize", scale.factor = 10000)


#identification of highly variable features (feature selection)
MUV27 <- FindVariableFeatures(MUV27, selection.method = "vst")


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(MUV27), 10)

# plot variable features with and without labels
#check to see if there are any immune related genes here
plot1 <- VariableFeaturePlot(MUV27)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
CombinePlots(plots = list(plot1, plot2))

#scaling the data
all.genes <- rownames(MUV27)
MUV27 <- ScaleData(MUV27, features = all.genes)

#perform PCA on scaled data using linear dimension reduction
#heat maps require PCA
MUV27 <- RunPCA(MUV27, features = VariableFeatures(object = MUV27), approx=FALSE)


#clustering cells
MUV27 <- FindNeighbors(MUV27, dims = 1:15)
MUV27 <- FindClusters(MUV27)

#determine marker groups and cluster data 
MUV27.markers <- FindAllMarkers(MUV27, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MUV27.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#find top markers and plot
top10 <- MUV27.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(MUV27) #+ theme(axis.text.y = element_text(size=3))

#want to extract the different clusters as individual matrices 
table(Idents(MUV27)) #how many cells are in each cluster 

# 0   1   2   --> cluster
# 111 110 93  --> number of cells per cluster

#extract cluster 0
cluster0.matrix.MUV27 <- as.matrix(GetAssayData(MUV27, slot = "counts")[, WhichCells(MUV27, ident = "0")])
ncol(cluster0.matrix.MUV27) #111 --> correct (column is cell)

#extract cluster 1
cluster1.matrix.MUV27 <- as.matrix(GetAssayData(MUV27, slot = "counts")[, WhichCells(MUV27, ident = "1")])
ncol(cluster1.matrix.MUV27) #110 --> correct  (column is cell)

#extract cluster 2
cluster2.matrix.MUV27 <- as.matrix(GetAssayData(MUV27, slot = "counts")[, WhichCells(MUV27, ident = "2")])
ncol(cluster2.matrix.MUV27) #93 --> correct  (column is cell)

#find the average expression value per cluster (sums all of the expression levels in the 
#cluster and divides by the number of cells --> average sum of expression per cell)

avg.c0.MUV27 <- sum(cluster0.matrix.MUV27)/ncol(cluster0.matrix.MUV27)
avg.c1.MUV27 <- sum(cluster1.matrix.MUV27)/ncol(cluster1.matrix.MUV27)
avg.c2.MUV27 <- sum(cluster2.matrix.MUV27)/ncol(cluster2.matrix.MUV27)

c("C0", "C1", "C2", "C3")[which.max(c(avg.c0.MUV27, avg.c1.MUV27, avg.c2.MUV27))] #C1


