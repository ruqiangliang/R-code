result <- NULL
therep <- function(name) {rep(name,5)}
cond <- lapply(c('saline24h','saline2w','ZNS24h','ZNS2w','control'),therep)
cond <- unlist(cond)
for (i in 1:24){
  result <- rbind(result, cbind(as.character(qPCR$Sample), as.character(qPCR$Gene), 
                                rep(cond[i],32),rep(colnames(qPCR)[i+2],32),qPCR[,i+2]))
}

##########ploting
result <- read.table('../Desktop/qPCR-result.txt',header=T)
qPCR.result <- summarySE(result, measurevar='expression',
                         groupvars=c('sample','gene','condition'))

## plot genes in organ of Corti
oc <- subset(qPCR.result, sample == 'OC')
oc <- subset(qPCR.result, sample == 'SGN')
genes <- unique(qPCR.result$gene)

i=1
  geneName <- as.character(genes[i])
  theGene <- subset(oc, gene == geneName)
  plot1 <- ggplot(theGene, aes(x=condition, y= expression)) +
    geom_errorbar(aes(ymin=expression-se, ymax=expression+se),width=.5)+
    geom_bar(stat='identity') +
    labs(title = geneName)

i=2
geneName <- as.character(genes[i])
theGene <- subset(oc, gene == geneName)
plot2 <- ggplot(theGene, aes(x=condition, y= expression)) +
  geom_errorbar(aes(ymin=expression-se, ymax=expression+se),width=.5)+
  geom_bar(stat='identity') +
  labs(title = geneName)

i=3
geneName <- as.character(genes[i])
theGene <- subset(oc, gene == geneName)
plot3 <- ggplot(theGene, aes(x=condition, y= expression)) +
  geom_errorbar(aes(ymin=expression-se, ymax=expression+se),width=.5)+
  geom_bar(stat='identity') +
  labs(title = geneName)

i=4
geneName <- as.character(genes[i])
theGene <- subset(oc, gene == geneName)
plot4 <- ggplot(theGene, aes(x=condition, y= expression)) +
  geom_errorbar(aes(ymin=expression-se, ymax=expression+se),width=.5)+
  geom_bar(stat='identity') +
  labs(title = geneName)

i=5
geneName <- as.character(genes[i])
theGene <- subset(oc, gene == geneName)
plot5 <- ggplot(theGene, aes(x=condition, y= expression)) +
  geom_errorbar(aes(ymin=expression-se, ymax=expression+se),width=.5)+
  geom_bar(stat='identity') +
  labs(title = geneName)

i=6
geneName <- as.character(genes[i])
theGene <- subset(oc, gene == geneName)
plot6 <- ggplot(theGene, aes(x=condition, y= expression)) +
  geom_errorbar(aes(ymin=expression-se, ymax=expression+se),width=.5)+
  geom_bar(stat='identity') +
  labs(title = geneName)

i=7
geneName <- as.character(genes[i])
theGene <- subset(oc, gene == geneName)
plot7 <- ggplot(theGene, aes(x=condition, y= expression)) +
  geom_errorbar(aes(ymin=expression-se, ymax=expression+se),width=.5)+
  geom_bar(stat='identity') +
  labs(title = geneName)

i=8
geneName <- as.character(genes[i])
theGene <- subset(oc, gene == geneName)
plot8 <- ggplot(theGene, aes(x=condition, y= expression)) +
  geom_errorbar(aes(ymin=expression-se, ymax=expression+se),width=.5)+
  geom_bar(stat='identity') +
  labs(title = geneName)

i=9
geneName <- as.character(genes[i])
theGene <- subset(oc, gene == geneName)
plot9 <- ggplot(theGene, aes(x=condition, y= expression)) +
  geom_errorbar(aes(ymin=expression-se, ymax=expression+se),width=.5)+
  geom_bar(stat='identity') +
  labs(title = geneName)

i=10
geneName <- as.character(genes[i])
theGene <- subset(oc, gene == geneName)
plot10 <- ggplot(theGene, aes(x=condition, y= expression)) +
  geom_errorbar(aes(ymin=expression-se, ymax=expression+se),width=.5)+
  geom_bar(stat='identity') +
  labs(title = geneName)

i=11
geneName <- as.character(genes[i])
theGene <- subset(oc, gene == geneName)
plot11 <- ggplot(theGene, aes(x=condition, y= expression)) +
  geom_errorbar(aes(ymin=expression-se, ymax=expression+se),width=.5)+
  geom_bar(stat='identity') +
  labs(title = geneName)

i=12
geneName <- as.character(genes[i])
theGene <- subset(oc, gene == geneName)
plot12 <- ggplot(theGene, aes(x=condition, y= expression)) +
  geom_errorbar(aes(ymin=expression-se, ymax=expression+se),width=.5)+
  geom_bar(stat='identity') +
  labs(title = geneName)

i=13
geneName <- as.character(genes[i])
theGene <- subset(oc, gene == geneName)
plot13 <- ggplot(theGene, aes(x=condition, y= expression)) +
  geom_errorbar(aes(ymin=expression-se, ymax=expression+se),width=.5)+
  geom_bar(stat='identity') +
  labs(title = geneName)

i=14
geneName <- as.character(genes[i])
theGene <- subset(oc, gene == geneName)
plot14 <- ggplot(theGene, aes(x=condition, y= expression)) +
  geom_errorbar(aes(ymin=expression-se, ymax=expression+se),width=.5)+
  geom_bar(stat='identity') +
  labs(title = geneName)

i=15
geneName <- as.character(genes[i])
theGene <- subset(oc, gene == geneName)
plot15 <- ggplot(theGene, aes(x=condition, y= expression)) +
  geom_errorbar(aes(ymin=expression-se, ymax=expression+se),width=.5)+
  geom_bar(stat='identity') +
  labs(title = geneName)

i=16
geneName <- as.character(genes[i])
theGene <- subset(oc, gene == geneName)
plot16 <- ggplot(theGene, aes(x=condition, y= expression)) +
  geom_errorbar(aes(ymin=expression-se, ymax=expression+se),width=.5)+
  geom_bar(stat='identity') +
  labs(title = geneName)

multiplot(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8,
          plot9, plot10, plot11, plot12, plot13, plot14, plot15, 
          plot16, cols=4)