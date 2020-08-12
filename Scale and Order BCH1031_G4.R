library(Seurat)
library("devtools")
library(Matrix)
library(dplyr)
library(scales)

#load in the data
BCH1031.data <- read.table("/Users/mayagarg/Downloads/GSE119926_RAW copy/GSM3905408_BCH1031.txt", header=TRUE)
#set up matrix
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

#scale each row 1 to 2 values 
for (i in 1:55) {
  x <- BCH1031mg.data[i,]!=0 #find and isolate nonzero values per row
  BCH1031mg.data[i,][x] <- scale(BCH1031mg.data[i,][x]) #scale nonzero values with each other 
  if(row.names(BCH1031mg.data[i,]) == "TMEM119") { BCH1031mg.data[i,][x] <- rescale(BCH1031mg.data[i,][x], to=c(2,3)) }
  else BCH1031mg.data[i,][x] <- rescale(BCH1031mg.data[i,][x], to=c(1,2)) # rescale down to preferred range 
}

#calculate the average of nonzero values per each cell and add to list
avg.vals.BCH1031 <- matrix(nrow = 273, ncol = 2)
for (y in 1:273) #loop through columns
  {nonzero <- 0
  sum <- 0 
  for (x in 1:55) #loop through rows
    {if (BCH1031mg.data[x,y] != 0)#if non zero
    { nonzero <- nonzero + 1  #increase nonzero count for the column
    sum <- sum + BCH1031mg.data[x,y]}}  #add value to sum for the column
  avg <- sum/nonzero#calculate the average for the column (I want to assign the average value to the cell name which is the column name)
  colnames(BCH1031mg.data[y]) <- avg#assign avg to coumn name 
  avg.vals.BCH1031[y, 1] <- (colnames(BCH1031mg.data[y])) #add that value to list
  avg.vals.BCH1031[y ,2] <- avg
}

avg.vals.BCH1031 <- avg.vals.BCH1031[order(avg.vals.BCH1031[,2], decreasing = TRUE), ]
avg.vals.BCH1031

