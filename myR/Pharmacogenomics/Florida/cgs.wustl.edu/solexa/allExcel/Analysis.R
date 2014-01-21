# load the libraries for CEPSKAT analysis
library('stringr')
library('CEPSKAT')

# compile the SNPs across all the samples together into one file
# befor that the Excel files were seperated into individual sheets by Gnumeric
setwd('/home/ruqiang/Pharacogenomics/Florida/cgs.wustl.edu/solexa/Bao_1062_2')
samples <- dir()[file.info(dir())$isdir]
for (i in 1:length(samples)){
  cmd <- paste('cp ', samples[i], '/snps_indels.xls ', '../allExcel/', 
               paste(samples[i],'_snps_indels.xls', sep=''),sep='')
  system(cmd)
}

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

# Read the SNP-Sample file and do CEPSKAT analysis

targetGT <- read.csv('Florida-CCB.csv',header=T)
# phenoT stores the phenotype information
phenoT <- read.csv('FloridaSampleP1-info.csv',header=T)
phenoT <- phenoT[phenoT$Barcode %in% unique(targetGT$ID),]

theP <- rep(NA, length(baitGenes))
theNum <- rep(0,length(baitGenes))

phenoT <- phenoT[order(phenoT$Z248),]
Y <- phenoT$Z248
nSample <- length(unique(targetGT$ID))
obj <- SKAT_Null_Model_CEP(Y ~ 1, Y[50], Y[49])

for (i in 1:length(baitGenes)) {
  res <- subset(targetGT, Gene == as.character(baitGenes[i]),select=c(1:5))
  if(length(res$ID) !=0){
    pos <- sort(unique(res$Position))
    matrixCol <- length(pos)
    theGT <- matrix(0,nSample, matrixCol)
    for (j in 1:nSample){
      theID <- phenoT$Barcode[j]
      snpGT <- subset(res,ID == as.character(theID),select=c(3:5))
      if (length(snpGT$Position) > 0) {
      thePos <- match(snpGT$Position, pos)
      for (m in 1:length(thePos)){
        theGT[j,thePos[m]] <- ifelse(str_detect(snpGT[m,3],toupper(snpGT[m,2])),1,2)
      }
    } }
    theP[i] <- SKAT(theGT, obj)$p.value
    theNum[i]<- matrixCol
    if (theP[i] < 0.05) {
      colnames(theGT)<-pos
      write.csv(cbind(phenoT,theGT),paste('./result/',baitGenes[i],'.csv',sep=''))
    }
  }
}
theReport <- cbind(as.character(baitGenes),theP, theNum)
colnames(theReport) <- c('Gene','P','Num')
write.csv(theReport,'./result/SKAT-report.csv', row.names=F)

geneCoverage <- NULL
for (i in 1:length(phenoT$Barcode)){
 theID <- phenoT$Barcode[i]
 res <- subset(targetGT, ID == as.character(theID), select=c('ID','Gene'))
 geneCoverage[i] <- length(unique(res$Gene))
}