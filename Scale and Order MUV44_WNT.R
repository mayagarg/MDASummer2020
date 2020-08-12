library(Seurat)
library("devtools")
library(Matrix)
library(dplyr)
library(scales)

#load in the data
MUV44.data <- read.table("/Users/mayagarg/Downloads/GSE119926_RAW copy/GSM3905418_MUV44.txt", header=TRUE)
#set up matrix
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

#scale each row 1 to 2 values 
for (i in 1:55) {
  x <- MUV44mg.data[i,]!=0 #find and isolate nonzero values per row
  MUV44mg.data[i,][x] <- scale(MUV44mg.data[i,][x]) #scale nonzero values with each other 
  if(row.names(MUV44mg.data[i,]) == "TMEM119") { MUV44mg.data[i,][x] <- rescale(MUV44mg.data[i,][x], to=c(2,3)) }
  else MUV44mg.data[i,][x] <- rescale(MUV44mg.data[i,][x], to=c(1,2)) # rescale down to preferred range 
}

#calculate the average of nonzero values per each cell and add to list
avg.vals.MUV44 <- matrix(nrow = 301, ncol = 2)
for (y in 1:301) #loop through columns
  {nonzero <- 0
  sum <- 0 
  for (x in 1:55) #loop through rows
    {if (MUV44mg.data[x,y] != 0)#if non zero
    { nonzero <- nonzero + 1  #increase nonzero count for the column
  sum <- sum + MUV44mg.data[x,y]}}  #add value to sum for the column
  avg <- sum/nonzero#calculate the average for the column (I want to assign the average value to the cell name which is the column name)
  colnames(MUV44mg.data[y]) <- avg#assign avg to coumn name 
  avg.vals.MUV44[y, 1] <- (colnames(MUV44mg.data[y])) #add that value to list
  avg.vals.MUV44[y ,2] <- avg
}

avg.vals.MUV44 <- avg.vals.MUV44[order(avg.vals.MUV44[,2], decreasing = TRUE), ]
avg.vals.MUV44

