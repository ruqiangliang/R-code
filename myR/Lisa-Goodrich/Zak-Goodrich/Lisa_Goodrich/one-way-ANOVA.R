myP <- NULL
myTemp <- tempfile(fileext = '.csv')
tissue <- c('SG_E12','SG_E13','SG_E16','SG_P0','SG_P6','SG_P15',
            + 'VG_E12','VG_E13','VG_E16','VG_P0','VG_P6','VG_P15')
theRep <- c(3,3,4,4,4,4,3,3,3,3,3,3)
myFrame <- rep(tissue,theRep)
for (i in 1:dim(gse)[1]) {
  myData <- gse[i,1:40]
  theFrame <- data.frame(cbind(myFrame, unlist(myData)))
  write.csv(theFrame,myTemp)
  theFrame <- read.csv(myTemp,header=T)
  theAOV <- aov(V2 ~ myFrame, data = theFrame)
  theANOVA <- anova(theAOV)
  myP[i] <- theANOVA[1,5]
}