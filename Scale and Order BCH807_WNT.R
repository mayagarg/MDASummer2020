library(Seurat)
library("devtools")
library(Matrix)
library(dplyr)
library(scales)

#load in the data
BCH807.data <- read.table("/Users/mayagarg/Downloads/GSE119926_RAW copy/GSM3905406_BCH807.txt", header=TRUE)
#set up matrix
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

#scale each row 1 to 2 values 
for (i in 1:55) {
  x <- BCH807mg.data[i,]!=0 #find and isolate nonzero values per row
  BCH807mg.data[i,][x] <- scale(BCH807mg.data[i,][x]) #scale nonzero values with each other 
  if(row.names(BCH807mg.data[i,]) == "TMEM119") { BCH807mg.data[i,][x] <- rescale(BCH807mg.data[i,][x], to=c(2,3)) }
  else BCH807mg.data[i,][x] <- rescale(BCH807mg.data[i,][x], to=c(1,2)) # rescale down to preferred range 
}

#calculate the average of nonzero values per each cell and add to list
avg.vals.BCH807 <- matrix(nrow = 307, ncol = 2)
for (y in 1:307) #loop through columns
  {nonzero <- 0
  sum <- 0 
  for (x in 1:55) #loop through rows
    {if (BCH807mg.data[x,y] != 0)#if non zero
    { nonzero <- nonzero + 1  #increase nonzero count for the column
    sum <- sum + BCH807mg.data[x,y]}}  #add value to sum for the column
  avg <- sum/nonzero#calculate the average for the column (I want to assign the average value to the cell name which is the column name)
  colnames(BCH807mg.data[y]) <- avg#assign avg to coumn name 
  avg.vals.BCH807[y, 1] <- (colnames(BCH807mg.data[y])) #add that value to list
  avg.vals.BCH807[y ,2] <- avg
}

avg.vals.BCH807 <- avg.vals.BCH807[order(avg.vals.BCH807[,2], decreasing = TRUE), ]
avg.vals.BCH807

