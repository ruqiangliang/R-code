myResult <- NULL
for (i in 1: length(gDNA$ID)){
  theID <- gDNA$ID[i]
  if (theID %in% allZ$V1){
  myResult <- rbind(myResult, cbind(gDNA[i,],subset(allZ,V1 == as.character(theID))))
}}