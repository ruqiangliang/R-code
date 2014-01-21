library('org.Mm.eg.db')

setwd("~/RNA-Seq/Analysis-R/")
allDir <- dir(pattern="Basal*")
for (m in 1: length(allDir)){
geneExp <- read.table(paste(allDir[m],'gene_exp.diff',sep='/'),header=T, stringsAsFactors=F)
gTable <- subset(geneExp, value_1 > 10 & value_2 >10 & abs(log2.fold_change.)>1 & gene != '-',
                 select=c('gene','value_1','value_2','log2.fold_change.'))
colnames(gTable) <- c('gene',paste('N',geneExp$sample_1[1]),paste('N',geneExp$sample_2[1]),'fold')
commaContain <- grep(',',gTable$gene)
for (i in 1: length(commaContain)){
  gnames <- strsplit(as.character(gTable$gene[commaContain[i]]),split=',')[[1]]
  gTable$gene[commaContain[i]] <- as.character(gnames[1])
  for (j in 2: length(gnames)){
    newRec <- cbind(gnames[j], gTable[commaContain[i],2:4])
    colnames(newRec)[1] <- 'gene'
    gTable <- rbind(gTable, newRec)
  }
}

x <- org.Mm.egSYMBOL2EG
# Get the entrez gene identifiers that are mapped to a gene symbol
mapped_genes <- mappedkeys(x)
# Convert to a list
symb2id <- as.list(x[mapped_genes])
theID <- symb2id[gTable$gene]

x <- org.Mm.egGENENAME
# Get the gene names that are mapped to an entrez gene identifier
mapped_genes <- mappedkeys(x)
# Convert to a list
id2name <- as.list(x[mapped_genes])
theNames <- NULL
for (i in 1: length(theID)){
  name <- id2name[unlist(theID[i])]
  theNames <- rbind(theNames, c(i,unlist(theID[i]), unlist(name)))
}
colnames(theNames) <- c('index','ID','Office.Name')
#theNames$id <- NULL
gTable <- cbind(gTable,theNames)
gTable <- gTable[order(gTable$fold),]
write.table(data.frame(gTable), paste(allDir[m],'-2013Dec09.txt',sep=''),sep='\t',row.names=F)
}