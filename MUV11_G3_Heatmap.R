library(Seurat)
library("devtools")
library(Matrix)
library(dplyr)

#load in the data
MUV11.data <- read.table("/Users/mayagarg/Downloads/GSE119926_RAW copy/GSM3905410_MUV11.txt", header=TRUE)

MUV11mg.data <- matrix() 
MUV11mg.data <- rbind(MUV11.data["SPARC", ], MUV11.data["P2RY12", ], MUV11.data["TMEM119", ], MUV11.data["SLC2A5", ], MUV11.data["SALL1", ], MUV11.data["ECSCR", ], 
                      MUV11.data["SCAMP5", ], MUV11.data["PTPRM", ], MUV11.data["KCND1", ], MUV11.data["SDK1", ], MUV11.data["CLSTN1", ], MUV11.data["SGCE", ], MUV11.data["SLC12A2", ],
                      MUV11.data["CADM1", ], MUV11.data["SEMA4G", ], MUV11.data["PROS1", ], MUV11.data["SLCO2B1", ], MUV11.data["RTN4RL1", ], MUV11.data["CMTM4", ], 
                      MUV11.data["FAT3", ], MUV11.data["SMO", ], MUV11.data["MRC2", ], MUV11.data["JAM2", ], MUV11.data["PLXNA4", ], MUV11.data["SLC46A1", ], 
                      MUV11.data["AGMO", ], MUV11.data["TMEM204", ], MUV11.data["BST2", ], MUV11.data["C2", ], MUV11.data["C4B", ], MUV11.data["CD74", ], 
                      MUV11.data["CFB", ], MUV11.data["CIITA", ], MUV11.data["CTSC", ], MUV11.data["EXO1", ], MUV11.data["GBP2", ], 
                      MUV11.data["IFNB1", ], MUV11.data["IRF1", ], MUV11.data["ITK", ], MUV11.data["JAK3", ], MUV11.data["LY86", ], MUV11.data["MASP1", ], 
                      MUV11.data["OAS2", ], MUV11.data["P2RY14", ], MUV11.data["SEMA4A", ], MUV11.data["SERPING1", ], MUV11.data["STAT1", ], 
                      MUV11.data["TLR2", ], MUV11.data["TLR3", ], MUV11.data["TNF", ], MUV11.data["TNFSF10", ], MUV11.data["TNFSF13b", ], MUV11.data["TNFSF15", ],MUV11.data["TNFSF8", ], 
                      MUV11.data["CX3CR1", ]) 


MUV11mg.data[is.na(MUV11mg.data)] <- 0

#initialize the Seurat data with the raw (non-normalized) data
#min.cells= include features (genes) detected in at least X many cells 
#min.features = include cells where at least Y many features (genes) are detected

MUV11 <- CreateSeuratObject(counts = MUV11mg.data, project= "GSM3905410_MUV11")
MUV11
#An object of class Seurat 
#55 features (genes) across 431 samples (cells) within 1 assay 
#Active assay: RNA (55 features, 0 variable features)

#QC pre processing
#MT is mitochondrial gene marker 
MUV11[["percent.mt"]] <- PercentageFeatureSet(MUV11, pattern = "^MT-")

#Normalizing data 
MUV11 <- NormalizeData(MUV11, normalization.method = "LogNormalize", scale.factor = 10000)


#identification of highly variable features (feature selection)
MUV11 <- FindVariableFeatures(MUV11, selection.method = "vst")

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(MUV11), 10)

# plot variable features with and without labels
#check to see if there are any immune related genes here
plot1 <- VariableFeaturePlot(MUV11)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
CombinePlots(plots = list(plot1, plot2))

#scaling the data
all.genes <- rownames(MUV11)
MUV11 <- ScaleData(MUV11, features = all.genes)

#perform PCA on scaled data using linear dimension reduction
#heat maps require PCA
MUV11 <- RunPCA(MUV11, features = VariableFeatures(object = MUV11), approx=FALSE)


#clustering cells
MUV11 <- FindNeighbors(MUV11, dims = 1:15)
MUV11 <- FindClusters(MUV11)

#determine marker groups and cluster data 
MUV11.markers <- FindAllMarkers(MUV11, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MUV11.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)



#find top markers and plot
top10 <- MUV11.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(MUV11) #+ theme(axis.text.y = element_text(size=3))

#want to extract the different clusters as individual matrices 
table(Idents(MUV11)) #how many cells are in each cluster 

# 0   1   2   --> cluster
# 174 134 123  --> number of cells per cluster

#extract cluster 0
cluster0.matrix.MUV11 <- as.matrix(GetAssayData(MUV11, slot = "counts")[, WhichCells(MUV11, ident = "0")])
ncol(cluster0.matrix.MUV11) #174 --> correct (column is cell)

#extract cluster 1
cluster1.matrix.MUV11 <- as.matrix(GetAssayData(MUV11, slot = "counts")[, WhichCells(MUV11, ident = "1")])
ncol(cluster1.matrix.MUV11) #134 --> correct  (column is cell)

#extract cluster 2
cluster2.matrix.MUV11 <- as.matrix(GetAssayData(MUV11, slot = "counts")[, WhichCells(MUV11, ident = "2")])
ncol(cluster2.matrix.MUV11) #123 --> correct  (column is cell)

#find the avergae expression value per cluster (sums all of the expression levels in the 
#cluster and divides by the number of cells --> average sum of expression per cell)

avg.c0.MUV11 <- sum(cluster0.matrix.MUV11)/ncol(cluster0.matrix.MUV11)
avg.c1.MUV11 <- sum(cluster1.matrix.MUV11)/ncol(cluster1.matrix.MUV11)
avg.c2.MUV11 <- sum(cluster2.matrix.MUV11)/ncol(cluster2.matrix.MUV11)

c("C0", "C1", "C2", "C3")[which.max(c(avg.c0.MUV11, avg.c1.MUV11, avg.c2.MUV11))] #C0


