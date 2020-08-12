library(Seurat)
library("devtools")
library(Matrix)
library(dplyr)
library(scales)

#load in the data
SJ454.data <- read.table("/Users/mayagarg/Downloads/GSE119926_RAW copy/GSM3905423_SJ454.txt", header=TRUE)
#set up matrix
SJ454mg.data <- matrix() 
SJ454mg.data <- rbind(SJ454.data["SPARC", ], SJ454.data["P2RY12", ], SJ454.data["TMEM119", ], SJ454.data["SLC2A5", ], SJ454.data["SALL1", ], SJ454.data["ECSCR", ], 
                      SJ454.data["SCAMP5", ], SJ454.data["PTPRM", ], SJ454.data["KCND1", ], SJ454.data["SDK1", ], SJ454.data["CLSTN1", ], SJ454.data["SGCE", ], SJ454.data["SLC12A2", ],
                      SJ454.data["CADM1", ], SJ454.data["SEMA4G", ], SJ454.data["PROS1", ], SJ454.data["SLCO2B1", ], SJ454.data["RTN4RL1", ], SJ454.data["CMTM4", ], 
                      SJ454.data["FAT3", ], SJ454.data["SMO", ], SJ454.data["MRC2", ], SJ454.data["JAM2", ], SJ454.data["PLXNA4", ], SJ454.data["SLC46A1", ], 
                      SJ454.data["AGMO", ], SJ454.data["TMEM204", ], SJ454.data["BST2", ], SJ454.data["C2", ], SJ454.data["C4B", ], SJ454.data["CD74", ], 
                      SJ454.data["CFB", ], SJ454.data["CIITA", ], SJ454.data["CTSC", ], SJ454.data["EXO1", ], SJ454.data["GBP2", ], 
                      SJ454.data["IFNB1", ], SJ454.data["IRF1", ], SJ454.data["ITK", ], SJ454.data["JAK3", ], SJ454.data["LY86", ], SJ454.data["MASP1", ], 
                      SJ454.data["OAS2", ], SJ454.data["P2RY14", ], SJ454.data["SEMA4A", ], SJ454.data["SERPING1", ], SJ454.data["STAT1", ], 
                      SJ454.data["TLR2", ], SJ454.data["TLR3", ], SJ454.data["TNF", ], SJ454.data["TNFSF10", ], SJ454.data["TNFSF13b", ], SJ454.data["TNFSF15", ],SJ454.data["TNFSF8", ], 
                      SJ454.data["CX3CR1", ]) 

SJ454mg.data[is.na(SJ454mg.data)] <- 0

#scale each row 1 to 2 values 
for (i in 1:55) {
  x <- SJ454mg.data[i,]!=0 #find and isolate nonzero values per row
  SJ454mg.data[i,][x] <- scale(SJ454mg.data[i,][x]) #scale nonzero values with each other 
  if(row.names(SJ454mg.data[i,]) == "TMEM119") { SJ454mg.data[i,][x] <- rescale(SJ454mg.data[i,][x], to=c(2,3)) }
  else SJ454mg.data[i,][x] <- rescale(SJ454mg.data[i,][x], to=c(1,2)) # rescale down to preferred range 
}

#calculate the average of nonzero values per each cell and add to list
avg.vals.SJ454 <- matrix(nrow = 293, ncol = 2)
for (y in 1:293) #loop through columns
  {nonzero <- 0
  sum <- 0 
  for (x in 1:55) #loop through rows
    {if (SJ454mg.data[x,y] != 0)#if non zero
    { nonzero <- nonzero + 1  #increase nonzero count for the column
    sum <- sum + SJ454mg.data[x,y]}}  #add value to sum for the column
  avg <- sum/nonzero#calculate the average for the column (I want to assign the average value to the cell name which is the column name)
  colnames(SJ454mg.data[y]) <- avg#assign avg to coumn name 
  avg.vals.SJ454[y, 1] <- (colnames(SJ454mg.data[y])) #add that value to list
  avg.vals.SJ454[y ,2] <- avg
}

avg.vals.SJ454 <- avg.vals.SJ454[order(avg.vals.SJ454[,2], decreasing = TRUE), ]
avg.vals.SJ454

