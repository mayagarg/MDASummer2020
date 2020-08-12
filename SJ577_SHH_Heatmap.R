library(Seurat)
library("devtools")
library(Matrix)
library(dplyr)

#load in the data
SJ577.data <- read.table("/Users/mayagarg/Downloads/GSE119926_RAW copy/GSM3905425_SJ577.txt", header=TRUE)

SJ577mg.data <- matrix() 
SJ577mg.data <- rbind(SJ577.data["SPARC", ], SJ577.data["P2RY12", ], SJ577.data["TMEM119", ], SJ577.data["SLC2A5", ], SJ577.data["SALL1", ], SJ577.data["ECSCR", ], 
                      SJ577.data["SCAMP5", ], SJ577.data["PTPRM", ], SJ577.data["KCND1", ], SJ577.data["SDK1", ], SJ577.data["CLSTN1", ], SJ577.data["SGCE", ], SJ577.data["SLC12A2", ],
                      SJ577.data["CADM1", ], SJ577.data["SEMA4G", ], SJ577.data["PROS1", ], SJ577.data["SLCO2B1", ], SJ577.data["RTN4RL1", ], SJ577.data["CMTM4", ], 
                      SJ577.data["FAT3", ], SJ577.data["SMO", ], SJ577.data["MRC2", ], SJ577.data["JAM2", ], SJ577.data["PLXNA4", ], SJ577.data["SLC46A1", ], 
                      SJ577.data["AGMO", ], SJ577.data["TMEM204", ], SJ577.data["BST2", ], SJ577.data["C2", ], SJ577.data["C4B", ], SJ577.data["CD74", ], 
                      SJ577.data["CFB", ], SJ577.data["CIITA", ], SJ577.data["CTSC", ], SJ577.data["EXO1", ], SJ577.data["GBP2", ], 
                      SJ577.data["IFNB1", ], SJ577.data["IRF1", ], SJ577.data["ITK", ], SJ577.data["JAK3", ], SJ577.data["LY86", ], SJ577.data["MASP1", ], 
                      SJ577.data["OAS2", ], SJ577.data["P2RY14", ], SJ577.data["SEMA4A", ], SJ577.data["SERPING1", ], SJ577.data["STAT1", ], 
                      SJ577.data["TLR2", ], SJ577.data["TLR3", ], SJ577.data["TNF", ], SJ577.data["TNFSF10", ], SJ577.data["TNFSF13b", ], SJ577.data["TNFSF15", ],SJ577.data["TNFSF8", ], 
                      SJ577.data["CX3CR1", ])

SJ577mg.data[is.na(SJ577mg.data)] <- 0

#initialize the Seurat data with the raw (non-normalized) data
#min.cells= include features (genes) detected in at least X many cells 
#min.features = include cells where at least Y many features (genes) are detected

SJ577 <- CreateSeuratObject(counts = SJ577mg.data, project= "GSM3905425_SJ577")
SJ577
#An object of class Seurat 
#55 features (genes) across 493 samples (cells) within 1 assay 
#Active assay: RNA (55 features, 0 variable features)

#QC pre processing
#MT is mitochondrial gene marker 
SJ577[["percent.mt"]] <- PercentageFeatureSet(SJ577, pattern = "^MT-")

#Normalizing data 
SJ577 <- NormalizeData(SJ577, normalization.method = "LogNormalize", scale.factor = 10000)


#identification of highly variable features (feature selection)
SJ577 <- FindVariableFeatures(SJ577, selection.method = "vst")


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(SJ577), 10)

# plot variable features with and without labels
#check to see if there are any immune related genes here
plot1 <- VariableFeaturePlot(SJ577)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
CombinePlots(plots = list(plot1, plot2))

#scaling the data
all.genes <- rownames(SJ577)
SJ577 <- ScaleData(SJ577, features = all.genes)

#perform PCA on scaled data using linear dimension reduction
#heat maps require PCA
SJ577 <- RunPCA(SJ577, features = VariableFeatures(object = SJ577), approx=FALSE)


#clustering cells
SJ577 <- FindNeighbors(SJ577, dims = 1:15)
SJ577 <- FindClusters(SJ577)

#determine marker groups and cluster data 
SJ577.markers <- FindAllMarkers(SJ577, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SJ577.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#find top markers and plot
top10 <- SJ577.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(SJ577) #+ theme(axis.text.y = element_text(size=3))

#want to extract the different clusters as individual matrices 
table(Idents(SJ577)) #how many cells are in each cluster 

# 0     1   2   3   4   --> cluster
# 138  131  87  75  62 --> number of cells per cluster

#extract cluster 0
cluster0.matrix.SJ577 <- as.matrix(GetAssayData(SJ577, slot = "counts")[, WhichCells(SJ577, ident = "0")])
ncol(cluster0.matrix.SJ577) #138 --> correct (column is cell)

#extract cluster 1
cluster1.matrix.SJ577 <- as.matrix(GetAssayData(SJ577, slot = "counts")[, WhichCells(SJ577, ident = "1")])
ncol(cluster1.matrix.SJ577) #131 --> correct  (column is cell)

#extract cluster 2
cluster2.matrix.SJ577 <- as.matrix(GetAssayData(SJ577, slot = "counts")[, WhichCells(SJ577, ident = "2")])
ncol(cluster2.matrix.SJ577) #87 --> correct  (column is cell)

#extract cluster 3
cluster3.matrix.SJ577 <- as.matrix(GetAssayData(SJ577, slot = "counts")[, WhichCells(SJ577, ident = "3")])
ncol(cluster3.matrix.SJ577) #75 --> correct  (column is cell)

#extract cluster 4
cluster4.matrix.SJ577 <- as.matrix(GetAssayData(SJ577, slot = "counts")[, WhichCells(SJ577, ident = "4")])
ncol(cluster4.matrix.SJ577) #62 --> correct  (column is cell)

avg.c0.SJ577 <- sum(cluster0.matrix.SJ577)/ncol(cluster0.matrix.SJ577)
avg.c1.SJ577 <- sum(cluster1.matrix.SJ577)/ncol(cluster1.matrix.SJ577)
avg.c2.SJ577 <- sum(cluster2.matrix.SJ577)/ncol(cluster2.matrix.SJ577)
avg.c3.SJ577 <- sum(cluster3.matrix.SJ577)/ncol(cluster3.matrix.SJ577)
avg.c4.SJ577 <- sum(cluster4.matrix.SJ577)/ncol(cluster4.matrix.SJ577)

c("C0", "C1", "C2", "C3", "C4")[which.max(c(avg.c0.SJ577, avg.c1.SJ577, avg.c2.SJ577, avg.c3.SJ577, avg.c4.SJ577))] #C4




