#library('utils')
library("pracma")
library("gdata")
library("plotrix")
library("VecStatGraphs3D")

rotationmat3D <- function (r,Axis) {
  L <- norm(Axis,"F")
  Axis <- Axis / L
  u <- Axis[1]; v <- Axis[2]; w <- Axis[3]
  u2 <- u^2; v2 <- v^2; w2 <- w^2
  c1 <- cos(r); s <- sin(r)
  R <- rep(NaN, times<-9); dim(R) <- c(3,3)
  R[1,1] <-  u2 + (v2 + w2)*c1; R[1,2] <- u*v*(1-c1) - w*s; R[1,3] <- u*w*(1-c1) + v*s
  R[2,1] <- u*v*(1-c1) + w*s; R[2,2] <- v2 + (u2+w2)*c1; R[2,3] <- v*w*(1-c1) - u*s
  R[3,1] <- u*w*(1-c1) - v*s; R[3,2] <- v*w*(1-c1)+u*s; R[3,3] <- w2+(u2+v2)*c1
  R
}
fname <- 'All2013Apr18_lowest'
#setwd(choose.dir())
setwd(tk_choose.dir(getwd(), "Choose a suitable folder"))
tf <- tempfile()
pdf(paste('../',fname,'Plot.pdf', sep='')) #All2012Nov29.pdf')
allData <- NULL
allIHC <- NULL

show(getwd())
noise <- dir('.')
fInfo <- file.info(noise)
dirIndex1 <- fInfo$isdir == TRUE
trueNoise <- noise[dirIndex1]

for (i in 1: sum(dirIndex1)) {
  setwd (file.path(getwd(), trueNoise[i], fsep= .Platform$file.sep))
  show(getwd())
  animal <- dir('.')
  fInfo <- file.info(animal)
  dirIndex2 <- fInfo$isdir == TRUE
  trueAnimal <- animal[dirIndex2]
  for (j in 1: sum(dirIndex2)) {
    setwd (file.path(getwd(), trueAnimal[j], fsep= .Platform$file.sep))
    show(getwd())
    freq <- dir('.')
    fInfo <- file.info(freq)
    dirIndex3 <- fInfo$isdir == TRUE
    trueFreq <- freq[dirIndex3]
    for (m in 1: sum(dirIndex3)) {
      setwd (file.path(getwd(), trueFreq[m], fsep= .Platform$file.sep))
      show(getwd())
      listCSV <- dir('.','*.csv')
      IHC2 <- read.csv(listCSV[2], stringsAsFactors=FALSE,header=FALSE)
      realLines <- c(IHC2$V8 == 'Cell' |  IHC2$V8 == 'Nuclei' |
                     IHC2$V8 == 'Ribbons' | IHC2$V8 == 'Receptors')
      IHC2 <- subset(IHC2[realLines,],select=c(paste('V', c(3,8,10,34:36),sep='')))
#      IHC2 <- data.frame(IHC2[realLines,],stringsAsFactors=F)
      colnames(IHC2) <- c('ID', 'Group', 'Volume','X','Y','Z')
#      IHC2$Volume <- as.numeric(IHC2$Volume)
#      IHC2$X <- as.numeric(IHC2$X)
#      IHC2$Y <- as.numeric(IHC2$Y)
#      IHC2$Z <- as.numeric(IHC2$Z)
#      synapses <- IHC2$Group == 'Ribbons' | IHC2$Group == 'Receptors'
#      lastBase <- c(IHC2$X[IHC2$Group == 'Nuclei'],
#                   IHC2$Y[IHC2$Group == 'Nuclei'],
#                   IHC2$Z[IHC2$Group == 'Nuclei'])
       lastBase <- subset(IHC2, Group == 'Nuclei',select = c('X','Y','Z'))
       lastBase <- as.numeric(lastBase)
      
      for (n in 1:length(listCSV)){
        ihc <- read.csv(listCSV[n], stringsAsFactors=FALSE,header=FALSE)
        realLines <- c( ihc$V8 == 'Cell' | ihc$V8 == 'Nuclei' |
              ihc$V8 == 'Receptors' | ihc$V8 == 'Ribbons')
        ihc <- subset(ihc[realLines,],select=c(paste('V', c(3,8,10,34:36),sep='')))
#        ihc <- data.frame(ihc)
        colnames(ihc) <- c('ID', 'Group', 'Volume','X','Y','Z')
#        ihc$Volume <- as.numeric(as.character(ihc$Volume))
#        ihc$X <- as.numeric(as.character(ihc$X))
#        ihc$Y <- as.numeric(as.character(ihc$Y))
#        ihc$Z <- as.numeric(as.character(ihc$Z))
        
        arrayRibbons <- ihc$Group == 'Ribbons'
        arrayReceptors <- ihc$Group == 'Receptors'
        synapses <- arrayRibbons | arrayReceptors
        
        if (n == 1) {
          theBase <- subset(ihc, Group == 'Nuclei', select = c('X','Y','Z'))
                       #c(ihc$X[ihc$Group == 'Nuclei'],
                       #ihc$Y[ihc$Group == 'Nuclei'],
                       #ihc$Z[ihc$Group == 'Nuclei'])
          xDirection <- 2*as.numeric(theBase) - as.numeric(lastBase)
        } else {
          xDirection <- as.numeric(theBase)
          theBase <- subset(ihc, Group == 'Nuclei', select = c('X','Y','Z')) #c(ihc$X[ihc$Group == 'Nuclei'],
                       #ihc$Y[ihc$Group == 'Nuclei'],
                       #ihc$Z[ihc$Group == 'Nuclei'])
        }
        dim(xDirection) <- c(1,3)
        theBase <- as.numeric(theBase)
        
        pointArray <- subset(ihc,select=c('X','Y','Z'))
        pointArray <- apply(pointArray,2,as.numeric)
        pointArray <- pointArray - rep(theBase,each=dim(pointArray)[1])
        
        write.table(pointArray,tf,row.names=F, col.names=F)
        dat <- LoadData3D(tf,Type=2)
        
        ribrec <- dat[synapses,1:6]
        radius <- ribrec[,1]
        quantRadius <- quantile(radius, probs = seq(0,1,0.2))        
        tips <- ribrec[radius > quantRadius[5],4:6]
        if (length(tips) >3){
          coordTip <- c(median(tips[,1]), median(tips[,2]), median(tips[,3]))}
        else {coordTip <- tips}
        
        #topCenter = c(ihc$X[ihc$Group == 'Cell'], ihc$Y[ihc$Group == 'Cell'],
        #              ihc$Z[ihc$Group == 'Cell'])
        #qx <- quantile(as.numeric(ihc$X[synapses]), probs = seq(0,1,0.1))
        #qy <- quantile(as.numeric(ihc$Y[synapses]), probs = seq(0,1,0.1))
        #qz <- quantile(as.numeric(ihc$Z[synapses]), probs = seq(0,1,0.1))
        #topCenter <- c((qx[2]+qx[10])/2, (qy[2]+qy[10])/2, (qz[2]+qz[10])/2)
        zDirection <- coordTip # topCenter - theBase
        dim(zDirection) <- c(1,3)
        if (zDirection[1,3] < 0){
          zDirection= -zDirection
        }
        
        
        zDirection <- zDirection/norm(zDirection,'F')
        theAxis <- matrix(c(0,0,1), nrow=1, ncol=3, byrow=TRUE)
        vtRotAxis <- cross(theAxis, zDirection)
        vtRotAxis <- vtRotAxis / norm(vtRotAxis, 'F')
        rotAngle <- acos(dot(theAxis, zDirection))
        tilt <- - rotAngle * 180/pi
        rotmat <- rotationmat3D(-rotAngle, vtRotAxis)
        pointArray <- pointArray %*% t(rotmat)
        arraySign <- sign(pointArray[,3])
        if (sum(arraySign == -1) < sum(arraySign == 1)) {
          pointArray <- - pointArray
        }
        
        
        vtXDir <- xDirection - theBase
        vtXDir <- vtXDir %*% t(rotmat) # vtXDir must be put into the new coordinate system
        rotAngle <- atan2(vtXDir[2],vtXDir[1])
        totateZ <- - rotAngle * 180/pi
        rotmat = rotationmat3D(-rotAngle, theAxis)
        pointArray <- pointArray %*% t(rotmat)
        #polArray <- cart2pol(pointArray)
        #polArray[,1] <- rad2deg(polArray[,1])
        
        
        write.table(pointArray,tf,row.names=F, col.names=F)
        dat <- LoadData3D(tf,Type=2)

        theAngle <- dat[,3]
        
        quad <- ifelse(theAngle < 45, 'Apical', ifelse(theAngle < 135,'Neural', 
                    ifelse(theAngle <225,'Basal', ifelse(theAngle < 315,'Abneural','Apical'))))
        quad2 <- ifelse(theAngle < 180, 'Neural2','Abneural2')
        RibbonSum <- c(sum(quad[arrayRibbons]=='Abneural'), sum(quad[arrayRibbons]=='Apical'),
                       sum(quad[arrayRibbons]=='Basal'), sum(quad[arrayRibbons]=='Neural'))
        RibbonSum2 <- c(sum(quad2[arrayRibbons]=='Abneural2'), sum(quad2[arrayRibbons]=='Neural2'))
        ReceptorSum <- c(sum(quad[arrayReceptors]=='Abneural'), sum(quad[arrayReceptors]=='Apical'),
                         sum(quad[arrayReceptors]=='Basal'), sum(quad[arrayReceptors]=='Neural'))
        ReceptorSum2 <- c(sum(quad2[arrayReceptors]=='Abneural2'), sum(quad2[arrayReceptors]=='Neural2'))
        allPoint <- dat[synapses,]
        thePoint <- allPoint[allPoint[,6]==min(allPoint[,6]),1:3]
        allIHC <- rbind(allIHC, c(trueNoise[i],trueAnimal[j],trueFreq[m],paste('IHC',n,sep=''),sum(RibbonSum),
                                  RibbonSum, RibbonSum / sum(RibbonSum),RibbonSum2, RibbonSum2/sum(RibbonSum),
                                  MeanDirection3D(dat[arrayRibbons,4:6]), MeanModule3D(dat[arrayRibbons,4:6]),
                                  sum(ReceptorSum),ReceptorSum,ReceptorSum/ sum(ReceptorSum), 
                                  ReceptorSum2,ReceptorSum2/ sum(ReceptorSum),
                                  MeanDirection3D(dat[arrayReceptors,4:6]),MeanModule3D(dat[arrayReceptors,4:6]), 
                                  thePoint))
        dd <- array()
        dd2 <- array()
        for (ii in 1: dim(dat[arrayRibbons,])[1]) {
          dd[ii] <- sort(sqrt((dat[arrayRibbons,4]-dat[arrayRibbons,4][ii])^2 + 
                           (dat[arrayRibbons,5]-dat[arrayRibbons,5][ii])^2+
                           (dat[arrayRibbons,6]-dat[arrayRibbons,6][ii])^2))[2]
        }
        
        for (ii in 1: dim(dat[arrayReceptors,])[1]) {
          dd2[ii] <- sort(sqrt((dat[arrayReceptors,4]-dat[arrayReceptors,4][ii])^2 + 
                                (dat[arrayReceptors,5]-dat[arrayReceptors,5][ii])^2+
                                (dat[arrayReceptors,6]-dat[arrayReceptors,6][ii])^2))[2]
        }

        ihc <- cbind(ihc, pointArray, dat[,1:3],quad,quad2,c(0,dd,dd2,0))
        infoIHC <- cbind(rep(trueNoise[i],times=length(ihc$X)), rep(trueAnimal[j],times=length(ihc$X)),
                         rep(trueFreq[m],times=length(ihc$X)), rep(paste('IHC',n,sep=''),times=length(ihc$X)))
        infoIHC <- cbind(infoIHC, ihc)
        allData <- rbind(allData,infoIHC)
        
        show(paste(trueNoise[i],trueAnimal[j],trueFreq[m],'IHC',n,sep='-'))
        
        ###########Plotting of the graph
        layout(matrix(c(1,1,2,3),2,2,byrow=T))
        polar.plot(dat[arrayRibbons,1],dat[arrayRibbons,3],rp.type='s', point.col='red',labels = c('Api','','Neur','','Bas','','Pil',''),
                   label.pos= c(0,45,90,135,180,225,270,315), radial.lim=c(0,max(dat[,1])),
                main=paste(paste(trueNoise[i],trueAnimal[j],trueFreq[m],'IHC',n,sep='-'),'X-Y, tilt: ', 
                format(tilt,digits=3), 'degree, rotate:', format(totateZ,digits=3),'degree'))
        polar.plot(dat[arrayReceptors,1],dat[arrayReceptors,3],rp.type='s',point.col='green',
                   radial.lim=c(0,max(polArray[,2])),add=T)
        
        polArray <- cart2pol(matrix(c(pointArray[,1], pointArray[,3]), nrow=length(ihc$X),ncol=2,byrow=F))
        polArray[,1] <- rad2deg(polArray[,1])
        polar.plot(polArray[arrayRibbons,2],polArray[arrayRibbons,1],rp.type='s', point.col='red', labels = c('Ap','','Top','','Bas','','Bottom',''),
                   label.pos= c(0,45,90,135,180,225,270,315),radial.lim=c(0,max(polArray[,2])),
                   main='X-Z')
        polar.plot(polArray[arrayReceptors,2],polArray[arrayReceptors,1],rp.type='s',point.col='green',
                   radial.lim=c(0,max(polArray[,2])),add=T)
        
        polArray <- cart2pol(matrix(c(pointArray[,2], pointArray[,3]), nrow=length(ihc$X),ncol=2,byrow=F))
        polArray[,1] <- rad2deg(polArray[,1])
        polar.plot(polArray[arrayRibbons,2],polArray[arrayRibbons,1],rp.type='s', point.col='red', labels = c('Neur','','Top','','Pil','','Bottom',''),
                   label.pos= c(0,45,90,135,180,225,270,315),radial.lim=c(0,max(polArray[,2])),
                   main='Y-Z')
        polar.plot(polArray[arrayReceptors,2],polArray[arrayReceptors,1],rp.type='s',point.col='green',
                   radial.lim=c(0,max(polArray[,2])),add=T)
      }
      
      setwd (file.path(getwd(), '..', fsep= .Platform$file.sep))  
    }
    setwd (file.path(getwd(), '..', fsep= .Platform$file.sep))
  }
  setwd (file.path(getwd(), '..', fsep= .Platform$file.sep))
}
#setwd (file.path(getwd(), '..', fsep= .Platform$file.sep))
allData <- data.frame(allData, stringsAsFactors=F)
allIHC <- data.frame(allIHC, stringsAsFactors=F)
colnames(allIHC) <- c('Noise','Animal', 'Frequency','IHC','Ribbon',
                      'RibbonAbneural','RibbonApical', 'RibbonBasal','RibbonNeural', 
                      'RibbonAbneuralPercent','RibbonApicalPercent','RibbonBasalPercent','RibbonNeuralPercent',
                      'RibbonAbneural2','RibbonNeural2','RibbonAbneural2Percent','RibbonNeural2Percent',
                      'RibbonColatitude','RibbonLongitude','RibbonModule','Receptor',
                      'ReceptorAbneural','ReceptorApical','ReceptorBasal','ReceptorNeural',
                      'ReceptorAbneuralPercent','ReceptorApicalPercent','ReceptorBasalPercent','ReceptorNeuralPercent',
                      'ReceptorAbneural2','ReceptorNeural2','ReceptorAbneural2Percent','ReceptorNeural2Percent',
                      'ReceptorColatitude','ReceptorLongitude','ReceptorModule',
                      'LowestMod','LowestColatitude','LowestLongitude')
colnames(allData) <- c('Noise','Animal', 'Frequency','IHC','ID','Group','Volume','X','Y','Z',
                       'CellX','CellY','CellZ','Module','Colatitude','Longitude','Quadrant','Quadrant2','MinDist')
write.csv(allData,paste('../',fname,'Synapse.csv', sep=''))
write.csv(allIHC,paste('../',fname,'Record.csv', sep=''))
dev.off()
unlink(tf)
#axis <- matrix(c(0,0,1),1,3, byrow=T)
#rotationmat3D(pi/4,axis)