library(Seurat)
library("devtools")
library(Matrix)
library(dplyr)
library(scales)

#load in the data
MUV29.data <- read.table("/Users/mayagarg/Downloads/GSE119926_RAW copy/GSM3905413_MUV29.txt", header=TRUE)
#set up matrix
MUV29mg.data <- matrix() 
MUV29mg.data <- rbind(MUV29.data["SPARC", ], MUV29.data["P2RY12", ], MUV29.data["TMEM119", ], MUV29.data["SLC2A5", ], MUV29.data["SALL1", ], MUV29.data["ECSCR", ], 
                       MUV29.data["SCAMP5", ], MUV29.data["PTPRM", ], MUV29.data["KCND1", ], MUV29.data["SDK1", ], MUV29.data["CLSTN1", ], MUV29.data["SGCE", ], MUV29.data["SLC12A2", ],
                       MUV29.data["CADM1", ], MUV29.data["SEMA4G", ], MUV29.data["PROS1", ], MUV29.data["SLCO2B1", ], MUV29.data["RTN4RL1", ], MUV29.data["CMTM4", ], 
                       MUV29.data["FAT3", ], MUV29.data["SMO", ], MUV29.data["MRC2", ], MUV29.data["JAM2", ], MUV29.data["PLXNA4", ], MUV29.data["SLC46A1", ], 
                       MUV29.data["AGMO", ], MUV29.data["TMEM204", ], MUV29.data["BST2", ], MUV29.data["C2", ], MUV29.data["C4B", ], MUV29.data["CD74", ], 
                       MUV29.data["CFB", ], MUV29.data["CIITA", ], MUV29.data["CTSC", ], MUV29.data["EXO1", ], MUV29.data["GBP2", ], 
                       MUV29.data["IFNB1", ], MUV29.data["IRF1", ], MUV29.data["ITK", ], MUV29.data["JAK3", ], MUV29.data["LY86", ], MUV29.data["MASP1", ], 
                       MUV29.data["OAS2", ], MUV29.data["P2RY14", ], MUV29.data["SEMA4A", ], MUV29.data["SERPING1", ], MUV29.data["STAT1", ], 
                       MUV29.data["TLR2", ], MUV29.data["TLR3", ], MUV29.data["TNF", ], MUV29.data["TNFSF10", ], MUV29.data["TNFSF13b", ], MUV29.data["TNFSF15", ],MUV29.data["TNFSF8", ], 
                       MUV29.data["CX3CR1", ]) 

MUV29mg.data[is.na(MUV29mg.data)] <- 0

#scale each row 1 to 2 values 
for (i in 1:55) {
  x <- MUV29mg.data[i,]!=0 #find and isolate nonzero values per row
  MUV29mg.data[i,][x] <- scale(MUV29mg.data[i,][x]) #scale nonzero values with each other 
  if(row.names(MUV29mg.data[i,]) == "TMEM119") { MUV29mg.data[i,][x] <- rescale(MUV29mg.data[i,][x], to=c(2,3)) }
  else MUV29mg.data[i,][x] <- rescale(MUV29mg.data[i,][x], to=c(1,2)) # rescale down to preferred range 
}

#calculate the average of nonzero values per each cell and add to list
avg.vals.MUV29 <- matrix(nrow = 314, ncol = 2)
for (y in 1:314) #loop through columns
  {nonzero <- 0
  sum <- 0 
    for (x in 1:55) #loop through rows
      {if (MUV29mg.data[x,y] != 0)#if non zero
      { nonzero <- nonzero + 1  #increase nonzero count for the column
      sum <- sum + MUV29mg.data[x,y]}}  #add value to sum for the column
    avg <- sum/nonzero#calculate the average for the column (I want to assign the average value to the cell name which is the column name)
    colnames(MUV29mg.data[y]) <- avg#assign avg to coumn name 
    avg.vals.MUV29[y, 1] <- (colnames(MUV29mg.data[y])) #add that value to list
    avg.vals.MUV29[y ,2] <- avg
}

avg.vals.MUV29 <- avg.vals.MUV29[order(avg.vals.MUV29[,2], decreasing = TRUE), ]
avg.vals.MUV29

