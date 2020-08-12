library(Seurat)
library("devtools")
library(Matrix)
library(dplyr)
library(scales)

#load in the data
MUV27.data <- read.table("/Users/mayagarg/Downloads/GSE119926_RAW copy/GSM3905412_MUV27.txt", header=TRUE)
#set up matrix
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

#scale each row 1 to 2 values 
for (i in 1:55) {
  x <- MUV27mg.data[i,]!=0 #find and isolate nonzero values per row
  MUV27mg.data[i,][x] <- scale(MUV27mg.data[i,][x]) #scale nonzero values with each other 
  if(row.names(MUV27mg.data[i,]) == "TMEM119") { MUV27mg.data[i,][x] <- rescale(MUV27mg.data[i,][x], to=c(2,3)) }
  else MUV27mg.data[i,][x] <- rescale(MUV27mg.data[i,][x], to=c(1,2)) # rescale down to preferred range 
}

#calculate the average of nonzero values per each cell and add to list
avg.vals.MUV27 <- matrix(nrow = 314, ncol = 2)
for (y in 1:314) #loop through columns
  {nonzero <- 0
  sum <- 0 
  for (x in 1:55) #loop through rows
    {if (MUV27mg.data[x,y] != 0)#if non zero
    { nonzero <- nonzero + 1  #increase nonzero count for the column
    sum <- sum + MUV27mg.data[x,y]}}  #add value to sum for the column
  avg <- sum/nonzero#calculate the average for the column (I want to assign the average value to the cell name which is the column name)
  colnames(MUV27mg.data[y]) <- avg#assign avg to coumn name 
  avg.vals.MUV27[y, 1] <- (colnames(MUV27mg.data[y])) #add that value to list
  avg.vals.MUV27[y ,2] <- avg
}

avg.vals.MUV27 <- avg.vals.MUV27[order(avg.vals.MUV27[,2], decreasing = TRUE), ]
avg.vals.MUV27

