############## Analysis of shuttle box data
# please e-mail the program, input and output
# to ruqiang.liang@gmail.com if you need me to change it
# No warranty is garrented.

library('stringr')
library('fields') # use image.plot function from fields package
allDATA <- NULL
myOS <- Sys.info()['sysname']  # choose.dir function is Windows specific
if (myOS == 'Linux' | myOS == 'Darwin'){
  require('tcltk')
  theDir <- tk_choose.dir(getwd(),caption='Please select the folder containing all the result file')
} else {
  theDir <- choose.dir(getwd(),'Please select the folder containing all the result file')
}

setwd(theDir)
pdf('../result-plot.pdf')
allFiles <- dir(theDir, pattern = '*.txt') # get all *.txt files from the shuttle box
mark <- 'TrainD'  # Filename is locked as XXXXXXXXTrainDxx.txt, before TrainD is mouseID, after is day
pos <- regexec(mark, allFiles) # extract mouseID and training day
for (i in 1:length(allFiles)){
  fname <- allFiles[i]  
  myd <- read.table(fname, comment='', header=T, sep='\t')
  numRec <- dim(myd)[1]
  mouseID <- substr(fname, 1, pos[[i]][1]-1)
  theDay <- substr(fname, pos[[i]][1]+str_length(mark), str_length(fname)-str_length('.txt'))
  ID <- rep(mouseID, numRec)
  Day <- rep(theDay, numRec)
  allDATA <- rbind(allDATA, cbind(ID, Day, myd))  # combine data from all the files in the folder into one dataframe
}
nd <- subset(allDATA, Trial. != 'Total:') # get rid of the Total: lines, only keeps the real data

totalRes <-NULL
idLevels <- unique(nd$ID) # seperate ID
for (id in 1:length(idLevels)){
  nnd <- subset(nd, ID == idLevels[id])
  dayLevels <- unique(nnd$Day)   # separate Day
  for (day in 1:length(dayLevels)){
    nnnd <- subset(nnd, Day == dayLevels[day])
    freqLevels <- unique(nnnd$kHz)   # seperate frequency
    for (freq in 1: length(freqLevels)){
      if (freqLevels[freq] != 0){   # frequency 0 is mostly passive
        nnnnd <- subset(nnnd, kHz == freqLevels[freq])
        totalActive <- sum(nnnnd$P.A == 'Active')
        rate <- sum(nnnnd$P.F == 'PASS')/totalActive   # calculate the pass rate for active test
        res <- c(as.character(idLevels[id]), as.character(dayLevels[day]), 
                 as.character(freqLevels[freq]), rate, totalActive)
        totalRes <- rbind(totalRes, res)   # combine all the results into one dataframe
      }
    }
  }
}
colnames(totalRes) <- c('ID','Day','kHz', 'Pass.Rate','Num.Active')
write.csv(totalRes,'../Result.csv',row.names=F)  # save as .csv file
totalRes <- read.csv('../Result.csv',header=T)   # read the csv file to convert into real numbers

for (id in 1:length(idLevels)){
  mouse <- subset(totalRes, ID == as.character(idLevels[id]))
  allFreq <- unique(mouse$kHz)
  allFreq <- sort(allFreq)    # sort the frequency for plot
  allDay <- unique(mouse$Day)
  allDay <- sort(allDay)      # sort the day for plot
  mouseMatrix <- matrix(NA, length(allDay),length(allFreq))    # construct a matrix for plot
  
  colnames(mouseMatrix)<- paste(allFreq,'kHz', sep='')         # set colnames as frequency
  rownames(mouseMatrix)<- paste('Day',allDay, sep='')          # set rownames as day
  totalMatrix <- mouseMatrix           # save another matrix for the total test number
  for (day in 1:length(allDay)){
    mouseDay <- subset(mouse, Day == allDay[day])
    mouseMatrix[day, match(mouseDay$kHz, allFreq)]<- mouseDay$Pass.Rate  # rewrite NA into real data
    totalMatrix[day, match(mouseDay$kHz, allFreq)]<- mouseDay$Num.Active # rewrite test number into real data
  }
  meanFreqMouse <- rowMeans(mouseMatrix, na.rm=T)
  freqMeanDay <- colMeans(mouseMatrix, na.rm=T)
  #### Plot the data
  image.plot(mouseMatrix, xlab='Day', ylab='Frequency',
             main=paste('Pass rate of ',as.character(idLevels[id]), sep=''))
  image.plot(totalMatrix, xlab='Day', ylab='Frequency',
                        main=paste('Active number of ',as.character(idLevels[id]), sep=''))
  plot(meanFreqMouse,type='b', xlab='Day', ylab='Active pass rate', 
       main=paste('Average pass rate of ',as.character(idLevels[id]), 
                  ' Across All Frequencies',sep=''), ylim=c(0,1))
  plot(freqMeanDay,type='b', xlab='Frequency (kHz)', ylab='Active pass rate', 
       main=paste('Average pass rate of ',as.character(idLevels[id]), 
                  ' Across All Days',sep=''), ylim=c(0,1), x=c(4:60))
  #### save the data
  mSaveMatrix <- rbind(cbind(mouseMatrix,meanFreqMouse), c(freqMeanDay,NA))
  write.csv(mSaveMatrix, paste('../Pass-Ratio-', as.character(idLevels[id]),'.csv', sep=''))
  write.csv(totalMatrix, paste('../Number-of-Active-', as.character(idLevels[id]),'.csv', sep=''))
}

## turn off the pdf printer
dev.off()
