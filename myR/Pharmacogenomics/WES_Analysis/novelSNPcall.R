targetGT <- NULL
snp <- dir('.','*.csv.2')
#baitGenes <- read.csv('./baitGenes.csv',header=T)
#baitGenes <- baitGenes$x
for (i in 1:length(snp)){
  snpF <- read.csv(snp[i],header=T)
  ID <- strsplit(strsplit(snp[i],split='_')[[1]][3],split='.csv.2')[[1]]
  GT <- subset(snpF, Gene %in% as.character(baitGenes))
  targetGT <- rbind(targetGT,c(ID, length(snpF$Gene),length(GT$Gene)))
}