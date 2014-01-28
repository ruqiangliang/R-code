# This script is to analyze the RNA-Sequencing data, four sample each group.
# The analysis is based on group data extracted from TopHat 2, cuffdiff function.
# It generates the cluster image
# By Ruqiang Liang, 2014

theData <- read.table('result.gct',skip=2,header=T)
newData <- subset(theData, (NoNoise != 0 | N90 !=0 | N96 != 0 | N96ZO != 0) & gene != '-')
theMatrix <- as.matrix(newData[,3:6])
colnames(theMatrix)<-c('No noise', '90 dB', '96 dB', '96dB ZNS')
rownames(theMatrix)<- newData$gene
theMean <- colMeans(theMatrix)
theFactor <- theMean/max(theMean)
newMatrix <- t(theMatrix)/theFactor # matrix normalized and transposed
theDist <- dist(newMatrix)
theHc <- hclust(theDist)
theDentr <- as.dendrogram(theHc)
pdf('Report-Bao-2014Jan27-cluster.pdf', paper='a4r')
par(mar=c(5,4,4,7)+0.1)
plot(theDentr,horiz=T,main='Cluster based on euclidean distance')
dev.off()
newData[,3:6] <- t(newMatrix)
NE90_over_NoNoise <- subset(newData, abs(log2(N90/NoNoise))>1 & (NoNoise * N90) != 0 & NoNoise > 10 & N90 > 10,
                            select = c('ID','gene','NoNoise','N90'))
NE90_over_NoNoise$log2_N90_over_NoNoise <- log2(NE90_over_NoNoise$N90/NE90_over_NoNoise$NoNoise)
NE90_over_NoNoise <- NE90_over_NoNoise[order(NE90_over_NoNoise$log2_N90_over_NoNoise),]

NE96_over_NoNoise <- subset(newData, abs(log2(N96/NoNoise))>1 & (NoNoise * N96) != 0 & NoNoise > 10 & N96 > 10,
                            select = c('ID','gene','NoNoise','N96'))
NE96_over_NoNoise$log2_N96_over_NoNoise <- log2(NE96_over_NoNoise$N96/NE96_over_NoNoise$NoNoise)
NE96_over_NoNoise <- NE96_over_NoNoise[order(NE96_over_NoNoise$log2_N96_over_NoNoise),]

NE96ZO_over_NoNoise <- subset(newData, abs(log2(N96ZO/NoNoise))>1 & (NoNoise * N96ZO) != 0 & NoNoise > 10 & N96ZO > 10,
                            select = c('ID','gene','NoNoise','N96ZO'))
NE96ZO_over_NoNoise$log2_N96ZO_over_NoNoise <- log2(NE96ZO_over_NoNoise$N96ZO/NE96ZO_over_NoNoise$NoNoise)
NE96ZO_over_NoNoise <- NE96ZO_over_NoNoise[order(NE96ZO_over_NoNoise$log2_N96ZO_over_NoNoise),]

NE96ZO_over_N96 <- subset(newData, abs(log2(N96ZO/N96))>1 & (N96 * N96ZO) != 0 & N96 > 10 & N96ZO > 10,
                              select = c('ID','gene','N96','N96ZO'))
NE96ZO_over_N96$log2_N96ZO_over_N96 <- log2(NE96ZO_over_N96$N96ZO/NE96ZO_over_N96$N96)
NE96ZO_over_N96 <- NE96ZO_over_N96[order(NE96ZO_over_N96$log2_N96ZO_over_N96),]

NE96ZO_over_N90 <- subset(newData, abs(log2(N96ZO/N90))>1 & (N90 * N96ZO) != 0 & N90 > 10 & N96ZO > 10,
                              select = c('ID','gene','N90','N96ZO'))
NE96ZO_over_N90$log2_N96ZO_over_N90 <- log2(NE96ZO_over_N90$N96ZO/NE96ZO_over_N90$N90)
NE96ZO_over_N90 <- NE96ZO_over_N90[order(NE96ZO_over_N90$log2_N96ZO_over_N90),]

NE96_over_N90 <- subset(newData, abs(log2(N96/N90))>1 & (N90 * N96) != 0 & N90 > 10 & N96 > 10,
                          select = c('ID','gene','N90','N96'))
NE96_over_N90$log2_N96_over_N90 <- log2(NE96_over_N90$N96/NE96_over_N90$N90)
NE96_over_N90 <- NE96_over_N90[order(NE96_over_N90$log2_N96_over_N90),]

#reports <- c(NE90_over_NoNoise, NE96_over_N90, NE96_over_NoNoise, NE96ZO_over_N90, 
#        NE96ZO_over_N96,NE96ZO_over_NoNoise)
write.table(NE90_over_NoNoise,'NE90_over_NoNoise.txt',row.names = F, sep='\t')
write.table(NE96_over_N90,'NE96_over_N90.txt',row.names = F, sep='\t')
write.table(NE96_over_NoNoise,'NE96_over_NoNoise.txt',row.names = F, sep='\t')
write.table(NE96ZO_over_N90,'NE96ZO_over_N90.txt',row.names = F, sep='\t')
write.table(NE96ZO_over_N96,'NE96ZO_over_N96.txt',row.names = F, sep='\t')
write.table(NE96ZO_over_NoNoise,'NE96ZO_over_NoNoise.txt',row.names = F, sep='\t')

NE96_high_over_NoNoise <- subset(NE96_over_NoNoise, log2_N96_over_NoNoise > 1)$gene
ZO96_down_over_NE96 <- subset(NE96ZO_over_N96, log2_N96ZO_over_N96 < -1)$gene
common_pressed_by_ZO <- intersect(NE96_high_over_NoNoise, ZO96_down_over_NE96)
common_pressed_by_ZO <- subset(newData, gene %in% common_pressed_by_ZO, select=c('ID','gene','NoNoise','N96','N96ZO'))
common_pressed_by_ZO$log2_N96_over_NoNoise <- log2(common_pressed_by_ZO$N96/common_pressed_by_ZO$NoNoise)
common_pressed_by_ZO$log2_N96ZO_over_N96 <- log2(common_pressed_by_ZO$N96ZO/common_pressed_by_ZO$N96)
write.table(common_pressed_by_ZO, 'Increase-in-N96-but-pressed-down-by-ZNS.txt', row.names=F, sep='\t')
real_pressed_by_ZO <- common_pressed_by_ZO

NE96_high_over_NoNoise <- subset(NE96_over_NoNoise, log2_N96_over_NoNoise < -1)$gene
ZO96_down_over_NE96 <- subset(NE96ZO_over_N96, log2_N96ZO_over_N96 > 1)$gene
common_pressed_by_ZO <- intersect(NE96_high_over_NoNoise, ZO96_down_over_NE96)
common_pressed_by_ZO <- subset(newData, gene %in% common_pressed_by_ZO, select=c('ID','gene','NoNoise','N96','N96ZO'))
common_pressed_by_ZO$log2_N96_over_NoNoise <- log2(common_pressed_by_ZO$N96/common_pressed_by_ZO$NoNoise)
common_pressed_by_ZO$log2_N96ZO_over_N96 <- log2(common_pressed_by_ZO$N96ZO/common_pressed_by_ZO$N96)
write.table(common_pressed_by_ZO, 'Decreased-in-N96-but-increase-by-ZNS.txt', row.names=F, sep='\t')
real_increased_by_ZO <- common_pressed_by_ZO
