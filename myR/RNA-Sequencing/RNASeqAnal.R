library(cummeRbund)
setwd("~/RNA-Seq")
cuff_data  <-readCufflinks('Basal_NoNoise_96dB')
csDensity(genes(cuff_data))
csScatter(genes(cuff_data), 'NoNoise', 'X96dB')
csVolcano(genes(cuff_data), 'NoNoise', 'X96dB')
mygene <-getGene(cuff_data,'Klf4')
expressionBarplot(mygene)
expressionBarplot(isoforms(mygene))

theResult <- NULL
for (i in 1:length(nonZeroTF$gene)){
  theDes <- subset(tf,Gene.Symbol == as.character(nonZeroTF$gene[i]))
  theResult <- rbind(theResult,theDes[1,])
}