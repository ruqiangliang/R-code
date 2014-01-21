
theDirs <- dir()
i=1
theData <- read.table(paste(theDirs[i],'/genes.fpkm_tracking', sep=''), header=T)
allData <- theData[,c(1,5,10,14)]
for (i in 2:length(theDirs)){
  theData <- read.table(paste(theDirs[i],'/genes.fpkm_tracking', sep=''), header=T)
  allData <- cbind(allData,theData[,c(10,14)])
}

write.table(allData, 'Combined-single-data.txt', row.names=F)
theColMean <- colMeans(allData[,3:18])
ScaleFactor <- theColMean / max(theColMean)
normData <- t(t(allData[,3:18]) / ScaleFactor)
normData <- cbind(allData[,1:2], normData)
write.table(normData , 'Combined-single-data-norm.txt', row.names=F)
afterData <- subset(normData, rowMeans(normData[,3:18]) != 0)
write.table(afterData, 'Combined-single-data-norm-rm0.gct', row.names=F, sep='\t')