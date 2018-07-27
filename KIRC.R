myData <- read.table("KIRC.txt", header = T, row.names = 1, stringsAsFactors = F)
myData <- as.matrix(myData)

#identifies the rownames for myData
rownames(myData)

#identifies dimensions of myData (row, col)
dim(myData)
# number of rows in myDatta
nrow(myData)
# number of Control columns
ncol(myData[,grep("Control", colnames(myData))])
# number of Tumor columns
ncol(myData[,grep("Tumor", colnames(myData))])

# creating log function from myData to use for the heatmap
myData_heatmap <- log((myData + 1), 2)
myData_heatmap

# saving the heatmap in a pdf
pdf(file = "Heatmap_KIRC_Data.pdf", width = 10, height = 10)
heatmap(myData_heatmap, 
         cluster_cols = F, 
         cluster_rows = T,
         scale = "row",
         clustering_distance_rows = "correlation",
         cellwidth = 25,
         cellheight = 10,
         color = colorRampPalette(c("blue", "white", "red"))(256),
         border_color = NA,
         main = "Heatmap of KIRC Data")
# ends the code in the pdf
dev.off()
dev.new()

#checks to see if myData is still numeric
is.numeric(myData)


tumorMean <- vector()
controlMean <- vector()
tumorsd <- vector()
controlsd <- vector()
pvals3 <- vector()
# dim(myData)[1] gets the number of columns in myData
for(k in 1:dim(myData)[1]){
  tumorMean <- c(tumorMean, mean(myData[k, grep("Tumor", colnames(myData))]))
  controlMean <- c(controlMean, mean(myData[k, grep("Control", colnames(myData))]))
  tumorsd <- c(tumorsd, sd(myData[k, grep("Tumor", colnames(myData))]))
  controlsd <- c(controlsd, sd(myData[k, grep("Control", colnames(myData))]))
  # right now it is a two-sided t-test, to make it one-sided, I have to add 
  # t.test(..., alternate = "greater")
  pval3 <- t.test(myData[k, grep("Tumor", colnames(myData))], 
                  myData[k, grep("Control", colnames(myData))])$p.value
  pvals3 <- c(pvals3, pval3)
}
# orders the tumor means by p-values smallest to greatest
tumorMean <- tumorMean[order(pvals3)]
# orders the control means by p-values smallest to greatest
controlMean <- controlMean[order(pvals3)]
# orders the tumor standard deviations by p-values smallest to greatest
tumorsd <- tumorsd[order(pvals3)]
# orders the control standard deviations by p-values smallest to greatest
controlsd <- controlsd[order(pvals3)]
# orders the miRNA by the p-values smallest to greatest
rowNames3 <- rownames(myData)[order(pvals3)]
# orders the p-values smallest to greatest
pvals3 <- pvals3[order(pvals3)]



# binds the vectors into columns
tempor <- cbind(tumorMean, controlMean, tumorsd, controlsd, pvals3)
# adds the miRNA as row names
row.names(tempor) = rowNames3
# makes sure tempor is a matrix
tempor <- as.matrix(tempor)
# checks that tempor is a matrix
is.matrix(tempor)
# gives the 30 most significant miRNAs
head(tempor, 30)

# createst a table of tempor data into a .txt
write.table(tempor, "KIRC_data.txt", sep = "\t", col.names = T)

#tempor= data.frame(tempor)
#is.data.frame(tempor)
#grep("138.9787754", tempor)
#tempor= data.frame(tempor)
#grep("138.97877538", tempor[,])

#write.table(cbind(tumorMean, controlMean, tumorsd, controlsd, pvals3),
            #"KIRC_data.txt", sep = "\t", col.names = colnames(("Tumor Mean",
            #"Control Mean", "Tumor SD", "Control SD", "P Value")))


