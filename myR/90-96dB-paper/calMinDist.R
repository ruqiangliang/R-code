distance <- read.csv('./Synapse-Analysis2012Nov29/AllRibbonReceptor.csv',sep='\t')
distRibbon <- subset(distance, Group=='Ribbons', select=c(Noise, Frequency, Animal,MinDist))
distReceptor <- subset(distance, Group=='Receptors', select=c(Noise, Frequency, Animal,MinDist))
noise <- unique (distance$Noise)
freq <- unique (distance$Frequency)
COM <- NULL

distance <- distRibbon
for (i in 1: length(noise)){
  for (j in 1: length(freq)){
    theData <- subset(distance, Noise == noise[i] & Frequency == freq[j], select = c(Animal,MinDist))
    animal <- unique(theData$Animal)
    for (m in 1: length(animal)){
      individual <- subset(theData, Animal == animal[m], select= MinDist)
      d <- density(individual$MinDist)
      Ind80 <- d$y > 0.2* max(d$y)
      theValue = with(d,sum(x[Ind80] * y[Ind80])/sum(y[Ind80]))
      COM <- rbind(COM,c(as.character(noise[i]), as.character(freq[j]), as.character(animal[m]),'Ribbons',
                         theValue))
    }
  }
}
distance <- distReceptor
for (i in 1: length(noise)){
  for (j in 1: length(freq)){
    theData <- subset(distance, Noise == noise[i] & Frequency == freq[j], select = c(Animal,MinDist))
    animal <- unique(theData$Animal)
    for (m in 1: length(animal)){
      individual <- subset(theData, Animal == animal[m], select= MinDist)
      d <- density(individual$MinDist)
      Ind80 <- d$y > 0.2* max(d$y)
      theValue = with(d,sum(x[Ind80] * y[Ind80])/sum(y[Ind80]))
      COM <- rbind(COM,c(as.character(noise[i]), as.character(freq[j]), as.character(animal[m]),'Receptors',
                         theValue))
    }
  }
}
colnames(COM) <- c('Noise','Freq','Animal','Group','MinDist')
rownames(COM) <- NULL