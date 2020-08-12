library(Seurat)
library("devtools")
library(Matrix)
library(dplyr)

#load in the data
MUV44.data <- read.table("/Users/mayagarg/Downloads/GSE119926_RAW copy/GSM3905418_MUV44.txt", header=TRUE)

MUV44mg.data <- matrix() 
MUV44mg.data <- rbind(MUV44.data["SPARC", ], MUV44.data["P2RY12", ], MUV44.data["TMEM119", ], MUV44.data["SLC2A5", ], MUV44.data["SALL1", ], MUV44.data["ECSCR", ], 
                      MUV44.data["SCAMP5", ], MUV44.data["PTPRM", ], MUV44.data["KCND1", ], MUV44.data["SDK1", ], MUV44.data["CLSTN1", ], MUV44.data["SGCE", ], MUV44.data["SLC12A2", ],
                      MUV44.data["CADM1", ], MUV44.data["SEMA4G", ], MUV44.data["PROS1", ], MUV44.data["SLCO2B1", ], MUV44.data["RTN4RL1", ], MUV44.data["CMTM4", ], 
                      MUV44.data["FAT3", ], MUV44.data["SMO", ], MUV44.data["MRC2", ], MUV44.data["JAM2", ], MUV44.data["PLXNA4", ], MUV44.data["SLC46A1", ], 
                      MUV44.data["AGMO", ], MUV44.data["TMEM204", ], MUV44.data["BST2", ], MUV44.data["C2", ], MUV44.data["C4B", ], MUV44.data["CD74", ], 
                      MUV44.data["CFB", ], MUV44.data["CIITA", ], MUV44.data["CTSC", ], MUV44.data["EXO1", ], MUV44.data["GBP2", ], 
                      MUV44.data["IFNB1", ], MUV44.data["IRF1", ], MUV44.data["ITK", ], MUV44.data["JAK3", ], MUV44.data["LY86", ], MUV44.data["MASP1", ], 
                      MUV44.data["OAS2", ], MUV44.data["P2RY14", ], MUV44.data["SEMA4A", ], MUV44.data["SERPING1", ], MUV44.data["STAT1", ], 
                      MUV44.data["TLR2", ], MUV44.data["TLR3", ], MUV44.data["TNF", ], MUV44.data["TNFSF10", ], MUV44.data["TNFSF13b", ], MUV44.data["TNFSF15", ],MUV44.data["TNFSF8", ], 
                      MUV44.data["CX3CR1", ])

MUV44mg.data[is.na(MUV44mg.data)] <- 0

#initialize the Seurat data with the raw (non-normalized) data
#min.cells= include features (genes) detected in at least X many cells 
#min.features = include cells where at least Y many features (genes) are detected

MUV44 <- CreateSeuratObject(counts = MUV44mg.data, project= "GSM3905418_MUV44")
MUV44
#An object of class Seurat 
#55 features (genes) across 301 samples (cells) within 1 assay 
#Active assay: RNA (55 features, 0 variable features)

#QC pre processing
#MT is mitochondrial gene marker 
MUV44[["percent.mt"]] <- PercentageFeatureSet(MUV44, pattern = "^MT-")

#Normalizing data 
MUV44 <- NormalizeData(MUV44, normalization.method = "LogNormalize", scale.factor = 10000)


#identification of highly variable features (feature selection)
MUV44 <- FindVariableFeatures(MUV44, selection.method = "vst")


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(MUV44), 10)

# plot variable features with and without labels
#check to see if there are any immune related genes here
plot1 <- VariableFeaturePlot(MUV44)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
CombinePlots(plots = list(plot1, plot2))

#scaling the data
all.genes <- rownames(MUV44)
MUV44 <- ScaleData(MUV44, features = all.genes)

#perform PCA on scaled data using linear dimension reduction
#heat maps require PCA
MUV44 <- RunPCA(MUV44, features = VariableFeatures(object = MUV44), approx=FALSE)

#clustering cells
MUV44 <- FindNeighbors(MUV44, dims = 1:15)
MUV44 <- FindClusters(MUV44)

#determine marker groups and cluster data 
MUV44.markers <- FindAllMarkers(MUV44, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MUV44.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#find top markers and plot
top10 <- MUV44.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(MUV44) #+ theme(axis.text.y = element_text(size=3))

#want to extract the different clusters as individual matrices 
table(Idents(MUV44)) #how many cells are in each cluster 

# 0   1   2   3 --> cluster
# 140 114 33  14  --> number of cells per cluster

#extract cluster 0
cluster0.matrix.MUV44 <- as.matrix(GetAssayData(MUV44, slot = "counts")[, WhichCells(MUV44, ident = "0")])
ncol(cluster0.matrix.MUV44) #140 --> correct (column is cell)

#extract cluster 1
cluster1.matrix.MUV44 <- as.matrix(GetAssayData(MUV44, slot = "counts")[, WhichCells(MUV44, ident = "1")])
ncol(cluster1.matrix.MUV44) #114 --> correct  (column is cell)

#extract cluster 2
cluster2.matrix.MUV44 <- as.matrix(GetAssayData(MUV44, slot = "counts")[, WhichCells(MUV44, ident = "2")])
ncol(cluster2.matrix.MUV44) #33 --> correct  (column is cell)

#extract cluster 3
cluster3.matrix.MUV44 <- as.matrix(GetAssayData(MUV44, slot = "counts")[, WhichCells(MUV44, ident = "3")])
ncol(cluster3.matrix.MUV44) #14 --> correct  (column is cell)

#find the average expression value per cluster (sums all of the expression levels in the 
#cluster and divides by the number of cells --> average sum of expression per cell)

avg.c0.MUV44 <- sum(cluster0.matrix.MUV44)/ncol(cluster0.matrix.MUV44)
avg.c1.MUV44 <- sum(cluster1.matrix.MUV44)/ncol(cluster1.matrix.MUV44)
avg.c2.MUV44 <- sum(cluster2.matrix.MUV44)/ncol(cluster2.matrix.MUV44)
avg.c3.MUV44 <- sum(cluster3.matrix.MUV44)/ncol(cluster3.matrix.MUV44)



c("C0", "C1", "C2", "C3")[which.max(c(avg.c0.MUV44, avg.c1.MUV44, avg.c2.MUV44, avg.c3.MUV44))] #C2


