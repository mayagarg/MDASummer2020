library(Seurat)
library("devtools")
library(Matrix)
library(dplyr)
library(scales)

#load in the data
SJ577.data <- read.table("/Users/mayagarg/Downloads/GSE119926_RAW copy/GSM3905425_SJ577.txt", header=TRUE)
#set up matrix
SJ577mg.data <- matrix() 
SJ577mg.data <- rbind(SJ577.data["SPARC", ], SJ577.data["P2RY12", ], SJ577.data["TMEM119", ], SJ577.data["SLC2A5", ], SJ577.data["SALL1", ], SJ577.data["ECSCR", ], 
                      SJ577.data["SCAMP5", ], SJ577.data["PTPRM", ], SJ577.data["KCND1", ], SJ577.data["SDK1", ], SJ577.data["CLSTN1", ], SJ577.data["SGCE", ], SJ577.data["SLC12A2", ],
                      SJ577.data["CADM1", ], SJ577.data["SEMA4G", ], SJ577.data["PROS1", ], SJ577.data["SLCO2B1", ], SJ577.data["RTN4RL1", ], SJ577.data["CMTM4", ], 
                      SJ577.data["FAT3", ], SJ577.data["SMO", ], SJ577.data["MRC2", ], SJ577.data["JAM2", ], SJ577.data["PLXNA4", ], SJ577.data["SLC46A1", ], 
                      SJ577.data["AGMO", ], SJ577.data["TMEM204", ], SJ577.data["BST2", ], SJ577.data["C2", ], SJ577.data["C4B", ], SJ577.data["CD74", ], 
                      SJ577.data["CFB", ], SJ577.data["CIITA", ], SJ577.data["CTSC", ], SJ577.data["EXO1", ], SJ577.data["GBP2", ], 
                      SJ577.data["IFNB1", ], SJ577.data["IRF1", ], SJ577.data["ITK", ], SJ577.data["JAK3", ], SJ577.data["LY86", ], SJ577.data["MASP1", ], 
                      SJ577.data["OAS2", ], SJ577.data["P2RY14", ], SJ577.data["SEMA4A", ], SJ577.data["SERPING1", ], SJ577.data["STAT1", ], 
                      SJ577.data["TLR2", ], SJ577.data["TLR3", ], SJ577.data["TNF", ], SJ577.data["TNFSF10", ], SJ577.data["TNFSF13b", ], SJ577.data["TNFSF15", ],SJ577.data["TNFSF8", ], 
                      SJ577.data["CX3CR1", ]) 

SJ577mg.data[is.na(SJ577mg.data)] <- 0

#scale each row 1 to 2 values 
for (i in 1:55) {
  x <- SJ577mg.data[i,]!=0 #find and isolate nonzero values per row
  SJ577mg.data[i,][x] <- scale(SJ577mg.data[i,][x]) #scale nonzero values with each other 
  if(row.names(SJ577mg.data[i,]) == "TMEM119") { SJ577mg.data[i,][x] <- rescale(SJ577mg.data[i,][x], to=c(2,3)) }
  else SJ577mg.data[i,][x] <- rescale(SJ577mg.data[i,][x], to=c(1,2)) # rescale down to preferred range 
}

#calculate the average of nonzero values per each cell and add to list
avg.vals.SJ577 <- matrix(nrow = 493, ncol = 2)
for (y in 1:493) #loop through columns
  {nonzero <- 0
  sum <- 0 
  for (x in 1:55) #loop through rows
    {if (SJ577mg.data[x,y] != 0)#if non zero
    { nonzero <- nonzero + 1  #increase nonzero count for the column
    sum <- sum + SJ577mg.data[x,y]}}  #add value to sum for the column
  avg <- sum/nonzero#calculate the average for the column (I want to assign the average value to the cell name which is the column name)
  colnames(SJ577mg.data[y]) <- avg#assign avg to coumn name 
  avg.vals.SJ577[y, 1] <- (colnames(SJ577mg.data[y])) #add that value to list
  avg.vals.SJ577[y ,2] <- avg
}

avg.vals.SJ577 <- avg.vals.SJ577[order(avg.vals.SJ577[,2], decreasing = TRUE), ]
avg.vals.SJ577

