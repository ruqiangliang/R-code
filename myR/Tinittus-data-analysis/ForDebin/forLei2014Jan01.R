########## Kernel density calculation of tinnitus data
# for any problem please e-mail the input, output, and the R program
# to ruqiang.liang@gmail.com
# (c) Ruqiang Liang

library(ks)

theDir <- choose.dir(getwd(), "Choose the folder containing the input csv files")
setwd(theDir)
allF <- dir('.', pattern = '*.csv')

for (fileF in 1: length(allF)){
  theF = allF[fileF]
fn <- basename(theF)

#pdf(paste(strsplit(fn,split='.csv'),'-plot.pdf',sep=''))

theRec <- NULL
raw <- read.csv(theF)
ID <- unique(raw$Encl) #levels(raw$EndlInfoA)
for (i in 1:length(ID)){
  mouse <- subset(raw, Encl == ID[i], select=c('TrialName','Max.N.'))
  trials <- levels(mouse$TrialName)
  freq <- c('4k', '8k','12.5k','16k','20k','25k')
  res <- NULL
  for (j in 1:length(freq)){
    gapPoint <- subset(mouse, TrialName == trials[grepl(freq[j], trials)& grepl ('ap', trials)])
    MaxN <- gapPoint$Max.N.
    theKDE <- kde(MaxN, hpi(MaxN),positive=T)
    Ind80 <- theKDE$estimate > 0.2*max(theKDE$estimate)
    gapValue <- sum(theKDE$estimate[Ind80]*theKDE$eval.points[Ind80])/sum(theKDE$estimate[Ind80])
    gapKDE <- theKDE
    gapMean <- mean(MaxN)
    
    startlePoint <- subset(mouse, TrialName == trials[grepl(freq[j], trials)& grepl ('tartle', trials)])
    MaxN <- startlePoint$Max.N.
    theKDE <- kde(MaxN, hpi(MaxN),positive=T)
    Ind80 <- theKDE$estimate > 0.2*max(theKDE$estimate)
    startleValue <- sum(theKDE$estimate[Ind80]*theKDE$eval.points[Ind80])/sum(theKDE$estimate[Ind80])
    startleKDE <- theKDE
    startleMean <- mean(MaxN)
    res <- c(res, gapValue, startleValue, gapValue/startleValue, gapMean, startleMean, gapMean/startleMean)
  }
  theRec <- cbind(theRec, res)
}

vlab <- c('gapKDE','startleKDE','ratioKDE', 'gapMean','startleMean','ratioMean')
myfun <- function(n,freq){
  f <- NULL
  for (i in 1: length(freq)){
    for (j in 1: length(n)){
      f <- c(f, paste(n[i],freq[j],sep='_'))
    }
  }
  f
}
rownames(theRec) <- myfun(freq,vlab)
colnames(theRec) <- seq(1:4)
  

write.csv(theRec,paste('../result/',strsplit(fn,split='.csv'),'-record.csv',sep=''))
}
#dev.off()