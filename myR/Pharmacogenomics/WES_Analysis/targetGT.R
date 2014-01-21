########### This program reads SNPs from sequencing result
# and phentoype data, then use CEPSKAT to calculate the association
# P value for each gene. If there is any problem, please e-mail the input,
# program, and output to Ruqiang Liang: ruqiang.liang@gmail.com
# (c) Ruqiang Liang, 2013. All rights reserved.


library('stringr')
library('CEPSKAT')

targetGT <- NULL
snp <- dir('.','*.csv.0')
baitGenes <- read.csv('./baitGenes.csv',header=T)
baitGenes <- unique(baitGenes$x)
for (i in 1:length(snp)){
  snpF <- read.csv(snp[i],header=T)
  ID <- strsplit(strsplit(snp[i],split='_')[[1]][1],split='.csv.0')[[1]]
  GT <- subset(snpF, Gene %in% as.character(baitGenes))
  targetGT <- rbind(targetGT, cbind(rep(ID,length(GT$Gene)),GT))
}

snp <- dir('.','*.csv.1')
#baitGenes <- read.csv('./baitGenes.csv',header=T)
#baitGenes <- unique(baitGenes$x)
for (i in 1:length(snp)){
  snpF <- read.csv(snp[i],header=T)
  ID <- strsplit(strsplit(snp[i],split='_')[[1]][1],split='.csv.1')[[1]]
  GT <- subset(snpF, Gene %in% as.character(baitGenes))
  targetGT <- rbind(targetGT, cbind(rep(ID,length(GT$Gene)),GT))
}

colnames(targetGT)[1] <- 'ID'
write.csv(targetGT,'Florida-CCB.csv',row.names=F)

phenoT <- read.csv('phenoT.csv',header=T)

theP <- rep(NA, length(baitGenes))
theNum <- rep(0,length(baitGenes))
phenoT <- phenoT[order(phenoT$Z),]
Y <- phenoT$Z

obj <- SKAT_Null_Model_CEP(Y ~ 1, Y[7], Y[6])

for (i in 1:length(baitGenes)){
  res <- subset(targetGT, Gene == as.character(baitGenes[i]),select=c(1:5))
  if(length(res$ID) !=0){
    pos <- sort(unique(res$Position))
    matrixCol <- length(pos)
    theGT <- matrix(0,12, matrixCol)
    for (j in 1:12){
      theID <- phenoT[j,2]
      snpGT <- subset(res,ID == theID,select=c(3:5))
      thePos <- match(snpGT$Position, pos)
      for (m in 1:length(thePos)){
        theGT[j,thePos[m]] <- ifelse(str_detect(snpGT[m,3],toupper(snpGT[m,2])),1,2)
      }
    }
    theP[i] <- SKAT(theGT, obj)$p.value
    theNum[i]<- matrixCol
    if (theP[i] < 0.25) {
      colnames(theGT)<-pos
      write.csv(cbind(phenoT,theGT),paste('./result2/',baitGenes[i],'.csv',sep=''))
    }
  }
}
theReport <- cbind(as.character(baitGenes),theP, theNum)
colnames(theReport) <- c('Gene','P','Num')
write.csv(theReport,'./result2/SKAT-report.csv')