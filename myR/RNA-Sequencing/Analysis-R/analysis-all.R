setwd("~/RNA-Seq/Analysis-R/")
allDir <- dir(pattern="Basal*")
allGeneExp <-list()
for (m in 1: length(allDir)){
  geneExp <- read.table(paste(allDir[m],'gene_exp.diff',sep='/'),header=T, stringsAsFactors=F)
  allGeneExp[[m]] <- geneExp
}
resultGeneExp <- list()
resultGeneExp$ID <- allGeneExp[[4]]$gene_id
resultGeneExp$gene <- allGeneExp[[4]]$gene
resultGeneExp$NoNoise <- allGeneExp[[4]]$value_1
resultGeneExp$N90 <- allGeneExp[[4]]$value_2
resultGeneExp$N96 <- allGeneExp[[5]]$value_2
resultGeneExp$N96ZO <- allGeneExp[[6]]$value_2
allResult <- as.data.frame(resultGeneExp)
write.table (allResult,'result.txt',row.names=F,sep='\t')