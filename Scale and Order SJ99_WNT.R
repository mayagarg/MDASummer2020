library(Seurat)
library("devtools")
library(Matrix)
library(dplyr)
library(scales)

#load in the data
SJ99.data <- read.table("/Users/mayagarg/Downloads/GSE119926_RAW copy/GSM3905420_SJ99.txt", header=TRUE)
#set up matrix
SJ99mg.data <- matrix() 
SJ99mg.data <- rbind(SJ99.data["SPARC", ], SJ99.data["P2RY12", ], SJ99.data["TMEM119", ], SJ99.data["SLC2A5", ], SJ99.data["SALL1", ], SJ99.data["ECSCR", ], 
                      SJ99.data["SCAMP5", ], SJ99.data["PTPRM", ], SJ99.data["KCND1", ], SJ99.data["SDK1", ], SJ99.data["CLSTN1", ], SJ99.data["SGCE", ], SJ99.data["SLC12A2", ],
                      SJ99.data["CADM1", ], SJ99.data["SEMA4G", ], SJ99.data["PROS1", ], SJ99.data["SLCO2B1", ], SJ99.data["RTN4RL1", ], SJ99.data["CMTM4", ], 
                      SJ99.data["FAT3", ], SJ99.data["SMO", ], SJ99.data["MRC2", ], SJ99.data["JAM2", ], SJ99.data["PLXNA4", ], SJ99.data["SLC46A1", ], 
                      SJ99.data["AGMO", ], SJ99.data["TMEM204", ], SJ99.data["BST2", ], SJ99.data["C2", ], SJ99.data["C4B", ], SJ99.data["CD74", ], 
                      SJ99.data["CFB", ], SJ99.data["CIITA", ], SJ99.data["CTSC", ], SJ99.data["EXO1", ], SJ99.data["GBP2", ], 
                      SJ99.data["IFNB1", ], SJ99.data["IRF1", ], SJ99.data["ITK", ], SJ99.data["JAK3", ], SJ99.data["LY86", ], SJ99.data["MASP1", ], 
                      SJ99.data["OAS2", ], SJ99.data["P2RY14", ], SJ99.data["SEMA4A", ], SJ99.data["SERPING1", ], SJ99.data["STAT1", ], 
                      SJ99.data["TLR2", ], SJ99.data["TLR3", ], SJ99.data["TNF", ], SJ99.data["TNFSF10", ], SJ99.data["TNFSF13b", ], SJ99.data["TNFSF15", ],SJ99.data["TNFSF8", ], 
                      SJ99.data["CX3CR1", ]) 

SJ99mg.data[is.na(SJ99mg.data)] <- 0

#scale each row 1 to 2 values 
for (i in 1:55) {
  x <- SJ99mg.data[i,]!=0 #find and isolate nonzero values per row
  SJ99mg.data[i,][x] <- scale(SJ99mg.data[i,][x]) #scale nonzero values with each other 
  if(row.names(SJ99mg.data[i,]) == "TMEM119") { SJ99mg.data[i,][x] <- rescale(SJ99mg.data[i,][x], to=c(2,3)) }
  else SJ99mg.data[i,][x] <- rescale(SJ99mg.data[i,][x], to=c(1,2)) # rescale down to preferred range 
}

#calculate the average of nonzero values per each cell and add to list
avg.vals.SJ99 <- matrix(nrow = 338, ncol = 2)
for (y in 1:339) #loop through columns
  {nonzero <- 0
  sum <- 0 
  for (x in 1:55) #loop through rows
    {if (SJ99mg.data[x,y] != 0)#if non zero
    { nonzero <- nonzero + 1  #increase nonzero count for the column
    sum <- sum + SJ99mg.data[x,y]}}  #add value to sum for the column
  avg <- sum/nonzero#calculate the average for the column (I want to assign the average value to the cell name which is the column name)
  colnames(SJ99mg.data[y]) <- avg#assign avg to coumn name 
  avg.vals.SJ99[y, 1] <- (colnames(SJ99mg.data[y])) #add that value to list
  avg.vals.SJ99[y ,2] <- avg
}

avg.vals.SJ99 <- avg.vals.SJ99[order(avg.vals.SJ99[,2], decreasing = TRUE), ]
avg.vals.SJ99

