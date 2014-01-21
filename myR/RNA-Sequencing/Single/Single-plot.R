#### Program to analysis individual RNA-Seq barcode
#### rather than in groups of barocoded samples
# This is calculating hierachyal cluster
# and draw the dendrograph
setwd("~/Desktop/BackUp/Projects/RNA-Sequencing/Single")
mydata <- read.table('Combined-single-data-norm-rm0.gct',skip=2, header=T)
mymatrix <- as.matrix(mydata[,3:18])
newname <- substr(oldname, 1, nchar(oldname)-5)
colnames(mymatrix) <- newname
mydist <- dist(t(mymatrix))
myhc <- hclust(mydist)
pdf('Single-hclust-param.pdf', paper='a4r')
par(mar=c(5,4,4,10)+0.1)
plot(as.dendrogram(myhc),horiz=T)
dev.off()
newdata <- subset(mydata, gene_short_name != '-')
mydata <- newdata
mymatrix <- as.matrix(mydata[,3:18])
newname <- substr(oldname, 1, nchar(oldname)-5)
colnames(mymatrix) <- newname
mydist <- dist(t(mymatrix))
myhc <- hclust(mydist)
pdf('Single-hclust-param-rmnull.pdf', paper='a4r')
par(mar=c(5,4,4,10)+0.1)
plot(as.dendrogram(myhc),horiz=T)
dev.off()