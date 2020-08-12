library(Seurat)
library("devtools")
library(Matrix)
library(dplyr)
library(scales)

#load in the data
BCH825.data <- read.table("/Users/mayagarg/Downloads/GSE119926_RAW copy/GSM3905407_BCH825.txt", header=TRUE)
#set up matrix
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

#scale each row 1 to 2 values 
for (i in 1:55) {
  x <- BCH825mg.data[i,]!=0 #find and isolate nonzero values per row
  BCH825mg.data[i,][x] <- scale(BCH825mg.data[i,][x]) #scale nonzero values with each other 
  if(row.names(BCH825mg.data[i,]) == "TMEM119") { BCH825mg.data[i,][x] <- rescale(BCH825mg.data[i,][x], to=c(2,3)) }
  else BCH825mg.data[i,][x] <- rescale(BCH825mg.data[i,][x], to=c(1,2)) # rescale down to preferred range 
}

#calculate the average of nonzero values per each cell and add to list
avg.vals.BCH825 <- matrix(nrow = 330, ncol = 2)
for (y in 1:330) #loop through columns
  {nonzero <- 0
  sum <- 0 
  for (x in 1:55) #loop through rows
    {if (BCH825mg.data[x,y] != 0)#if non zero
    { nonzero <- nonzero + 1  #increase nonzero count for the column
      sum <- sum + BCH825mg.data[x,y]}}  #add value to sum for the column
  avg <- sum/nonzero#calculate the average for the column (I want to assign the average value to the cell name which is the column name)
  colnames(BCH825mg.data[y]) <- avg#assign avg to coumn name 
  avg.vals.BCH825[y, 1] <- (colnames(BCH825mg.data[y])) #add that value to list
  avg.vals.BCH825[y ,2] <- avg
  }

avg.vals.BCH825 <- avg.vals.BCH825[order(avg.vals.BCH825[,2], decreasing = TRUE), ]
avg.vals.BCH825

