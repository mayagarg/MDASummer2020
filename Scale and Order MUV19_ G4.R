library(Seurat)
library("devtools")
library(Matrix)
library(dplyr)
library(scales)

#load in the data
MUV19.data <- read.table("/Users/mayagarg/Downloads/GSE119926_RAW copy/GSM3905411_MUV19.txt", header=TRUE)
#set up matrix
MUV19mg.data <- matrix() 
MUV19mg.data <- rbind(MUV19.data["SPARC", ], MUV19.data["P2RY12", ], MUV19.data["TMEM119", ], MUV19.data["SLC2A5", ], MUV19.data["SALL1", ], MUV19.data["ECSCR", ], 
                      MUV19.data["SCAMP5", ], MUV19.data["PTPRM", ], MUV19.data["KCND1", ], MUV19.data["SDK1", ], MUV19.data["CLSTN1", ], MUV19.data["SGCE", ], MUV19.data["SLC12A2", ],
                      MUV19.data["CADM1", ], MUV19.data["SEMA4G", ], MUV19.data["PROS1", ], MUV19.data["SLCO2B1", ], MUV19.data["RTN4RL1", ], MUV19.data["CMTM4", ], 
                      MUV19.data["FAT3", ], MUV19.data["SMO", ], MUV19.data["MRC2", ], MUV19.data["JAM2", ], MUV19.data["PLXNA4", ], MUV19.data["SLC46A1", ], 
                      MUV19.data["AGMO", ], MUV19.data["TMEM204", ], MUV19.data["BST2", ], MUV19.data["C2", ], MUV19.data["C4B", ], MUV19.data["CD74", ], 
                      MUV19.data["CFB", ], MUV19.data["CIITA", ], MUV19.data["CTSC", ], MUV19.data["EXO1", ], MUV19.data["GBP2", ], 
                      MUV19.data["IFNB1", ], MUV19.data["IRF1", ], MUV19.data["ITK", ], MUV19.data["JAK3", ], MUV19.data["LY86", ], MUV19.data["MASP1", ], 
                      MUV19.data["OAS2", ], MUV19.data["P2RY14", ], MUV19.data["SEMA4A", ], MUV19.data["SERPING1", ], MUV19.data["STAT1", ], 
                      MUV19.data["TLR2", ], MUV19.data["TLR3", ], MUV19.data["TNF", ], MUV19.data["TNFSF10", ], MUV19.data["TNFSF13b", ], MUV19.data["TNFSF15", ],MUV19.data["TNFSF8", ], 
                      MUV19.data["CX3CR1", ]) 

MUV19mg.data[is.na(MUV19mg.data)] <- 0

#scale each row 1 to 2 values 
for (i in 1:55) {
  x <- MUV19mg.data[i,]!=0 #find and isolate nonzero values per row
  MUV19mg.data[i,][x] <- scale(MUV19mg.data[i,][x]) #scale nonzero values with each other 
  if(row.names(MUV19mg.data[i,]) == "TMEM119") { MUV19mg.data[i,][x] <- rescale(MUV19mg.data[i,][x], to=c(2,3)) }
  else MUV19mg.data[i,][x] <- rescale(MUV19mg.data[i,][x], to=c(1,2)) # rescale down to preferred range 
}

#calculate the average of nonzero values per each cell and add to list
avg.vals.MUV19 <- matrix(nrow = 304, ncol = 2)
for (y in 1:304) #loop through columns
  {nonzero <- 0
  sum <- 0 
  for (x in 1:55) #loop through rows
    {if (MUV19mg.data[x,y] != 0)#if non zero
    { nonzero <- nonzero + 1  #increase nonzero count for the column
    sum <- sum + MUV19mg.data[x,y]}}  #add value to sum for the column
  avg <- sum/nonzero#calculate the average for the column (I want to assign the average value to the cell name which is the column name)
  colnames(MUV19mg.data[y]) <- avg#assign avg to coumn name 
  avg.vals.MUV19[y, 1] <- (colnames(MUV19mg.data[y])) #add that value to list
  avg.vals.MUV19[y ,2] <- avg
}

avg.vals.MUV19 <- avg.vals.MUV19[order(avg.vals.MUV19[,2], decreasing = TRUE), ]
avg.vals.MUV19

