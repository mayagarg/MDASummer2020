library(Seurat)
library("devtools")
library(Matrix)
library(dplyr)
library(scales)

#load in the data
MUV41.data <- read.table("/Users/mayagarg/Downloads/GSE119926_RAW copy/GSM3905417_MUV41.txt", header=TRUE)
#set up matrix
MUV41mg.data <- matrix() 
MUV41mg.data <- rbind(MUV41.data["SPARC", ], MUV41.data["P2RY12", ], MUV41.data["TMEM119", ], MUV41.data["SLC2A5", ], MUV41.data["SALL1", ], MUV41.data["ECSCR", ], 
                      MUV41.data["SCAMP5", ], MUV41.data["PTPRM", ], MUV41.data["KCND1", ], MUV41.data["SDK1", ], MUV41.data["CLSTN1", ], MUV41.data["SGCE", ], MUV41.data["SLC12A2", ],
                      MUV41.data["CADM1", ], MUV41.data["SEMA4G", ], MUV41.data["PROS1", ], MUV41.data["SLCO2B1", ], MUV41.data["RTN4RL1", ], MUV41.data["CMTM4", ], 
                      MUV41.data["FAT3", ], MUV41.data["SMO", ], MUV41.data["MRC2", ], MUV41.data["JAM2", ], MUV41.data["PLXNA4", ], MUV41.data["SLC46A1", ], 
                      MUV41.data["AGMO", ], MUV41.data["TMEM204", ], MUV41.data["BST2", ], MUV41.data["C2", ], MUV41.data["C4B", ], MUV41.data["CD74", ], 
                      MUV41.data["CFB", ], MUV41.data["CIITA", ], MUV41.data["CTSC", ], MUV41.data["EXO1", ], MUV41.data["GBP2", ], 
                      MUV41.data["IFNB1", ], MUV41.data["IRF1", ], MUV41.data["ITK", ], MUV41.data["JAK3", ], MUV41.data["LY86", ], MUV41.data["MASP1", ], 
                      MUV41.data["OAS2", ], MUV41.data["P2RY14", ], MUV41.data["SEMA4A", ], MUV41.data["SERPING1", ], MUV41.data["STAT1", ], 
                      MUV41.data["TLR2", ], MUV41.data["TLR3", ], MUV41.data["TNF", ], MUV41.data["TNFSF10", ], MUV41.data["TNFSF13b", ], MUV41.data["TNFSF15", ],MUV41.data["TNFSF8", ], 
                      MUV41.data["CX3CR1", ]) 

MUV41mg.data[is.na(MUV41mg.data)] <- 0

#scale each row 1 to 2 values 
for (i in 1:55) {
  x <- MUV41mg.data[i,]!=0 #find and isolate nonzero values per row
  MUV41mg.data[i,][x] <- scale(MUV41mg.data[i,][x]) #scale nonzero values with each other 
  if(row.names(MUV41mg.data[i,]) == "TMEM119") { MUV41mg.data[i,][x] <- rescale(MUV41mg.data[i,][x], to=c(2,3)) }
  else MUV41mg.data[i,][x] <- rescale(MUV41mg.data[i,][x], to=c(1,2)) # rescale down to preferred range 
}

#calculate the average of nonzero values per each cell and add to list
avg.vals.MUV41 <- matrix(nrow = 338, ncol = 2)
for (y in 1:338) #loop through columns
  {nonzero <- 0
  sum <- 0 
  for (x in 1:55) #loop through rows
    {if (MUV41mg.data[x,y] != 0)#if non zero
    { nonzero <- nonzero + 1  #increase nonzero count for the column
  sum <- sum + MUV41mg.data[x,y]}}  #add value to sum for the column
  avg <- sum/nonzero#calculate the average for the column (I want to assign the average value to the cell name which is the column name)
  colnames(MUV41mg.data[y]) <- avg#assign avg to coumn name 
  avg.vals.MUV41[y, 1] <- (colnames(MUV41mg.data[y])) #add that value to list
  avg.vals.MUV41[y ,2] <- avg
}

avg.vals.MUV41 <- avg.vals.MUV41[order(avg.vals.MUV41[,2], decreasing = TRUE), ]
avg.vals.MUV41

