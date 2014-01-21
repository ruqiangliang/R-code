library('stringr')
library('AssotesteR')
theP <- list()
### Test list is required by the MULTI function
# Pay attention to RVT1 method. It was used by the cystic fibrosis paper in Nature Genetics 44(2012):886
test_list = c("BST", "CMAT", "CALPHA", "ORWSS", "RWAS", "RBT", "SCORE",
              "SUM", "SSU", "WSS", "WST","RVT1","RVT2","VT")
phenoT <- phenoT[order(phenoT$Z248),]
Y <- phenoT$Z248
nSample <- length(unique(targetGT$ID))
# The multi association test calculate on dichotimous phenotype rather than continuous extreme phenotype
# so we need to cast continous phentoype into dichotimous phenotype
newY <- Y
newY[newY < 0] <- 0
newY[newY > 0] <- 1
# again, construct the genotype matrix for each gene
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
    if (dim(theGT)[2]>1){
      # We use the MULTI function from the AssotesteR package to perform the association test
      # We set the minor allele frequency cut off value 0.125 which was used in Nature Genetics 44(2012):886
      # And do 500 permutation
      mymulti1 <- MULTI(newY, theGT, test_list, maf=0.125, perm= 500)     
      # Save the pvalue of the association test methods listed in the list
      theRow <- c(as.character(baitGenes[i]),matrixCol, mymulti1$pvalue)
      theP[[i]] <- theRow  # save the result into a list
      show(theRow) # show the result to the console so that you can see what's going on
    }
  }
}
# Change list into dataframe so we can save into a file
theReport <- NULL
for (i in 1:length(theP)){
  theReport <- rbind(theReport, theP[[i]])
}

colnames(theReport) <- c('Gene','Num', row.names(mymulti1))
write.csv(theReport,'MULTI-report.csv', row.names=F)