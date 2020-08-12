library(Seurat)
library("devtools")
library(Matrix)
library(dplyr)
library(scales)

#load in the data
MUV11.data <- read.table("/Users/mayagarg/Downloads/GSE119926_RAW copy/GSM3905410_MUV11.txt", header=TRUE)
#set up matrix
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

#scale each row 1 to 2 values 
for (i in 1:55) {
  x <- MUV11mg.data[i,]!=0 #find and isolate nonzero values per row
  MUV11mg.data[i,][x] <- scale(MUV11mg.data[i,][x]) #scale nonzero values with each other 
  if(row.names(MUV11mg.data[i,]) == "TMEM119") { MUV11mg.data[i,][x] <- rescale(MUV11mg.data[i,][x], to=c(2,3)) }
  else MUV11mg.data[i,][x] <- rescale(MUV11mg.data[i,][x], to=c(1,2)) # rescale down to preferred range 
}

#calculate the average of nonzero values per each cell and add to list
avg.vals.MUV11 <- matrix(nrow = 431, ncol = 2)
for (y in 1:431) #loop through columns
  {nonzero <- 0
  sum <- 0 
for (x in 1:55) #loop through rows
  {if (MUV11mg.data[x,y] != 0)#if non zero
    { nonzero <- nonzero + 1  #increase nonzero count for the column
      sum <- sum + MUV11mg.data[x,y]}}  #add value to sum for the column
  avg <- sum/nonzero#calculate the average for the column (I want to assign the average value to the cell name which is the column name)
  colnames(MUV11mg.data[y]) <- avg #assign avg to coumn name 
  avg.vals.MUV11[y, 1] <- (colnames(MUV11mg.data[y])) #add that value to list
  avg.vals.MUV11[y ,2] <- avg
}

avg.vals.MUV11 <- avg.vals.MUV11[order(avg.vals.MUV11[,2], decreasing = TRUE), ]
avg.vals.MUV11

