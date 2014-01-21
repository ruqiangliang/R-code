
  pos <- orderSNP
  matrixCol <- length(pos)
  theGT <- matrix(0,12, matrixCol*2)
  for (i in 1:matrixCol){
    theGT[,(i*2-1):(i*2)] <- rep(toupper(myMap[i,1]),2)
  }
  
  for (j in 1:12){
    theID <- phenoT[j,2]
    snpGT <- subset(temp,ID == theID,select=c(5,7))
    thePos <- match(snpGT$allSNP, orderSNP)
    for (m in 1:length(thePos)){
      theGT[j,(2*thePos[m]-1):(2*thePos[m])] <- strsplit(as.character(snpGT[m,1])
                                                         ,split='/')[[1]]
    }
  }

  myMap <- NULL
  for (i in 1: length(orderSNP)){
    snp <- orderSNP[i]
    chr <- strsplit(snp,split='p')[[1]]
    myMap <- rbind(myMap,c(str_sub(chr[1],4), snp, 0, chr[2]))
  }
  
  myMap <- NULL
  for (i in 1: length(orderSNP)){
    snp <- orderSNP[i]
    myMap <- rbind(myMap,subset(tempSet, allSNP == as.character(orderSNP[i]))[1,])
  }