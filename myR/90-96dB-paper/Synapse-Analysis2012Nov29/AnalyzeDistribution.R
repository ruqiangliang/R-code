require(ggplot2)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE, conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This is does the summary; it's not easy to understand...
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun= function(xx, col, na.rm) {
                   c( N    = length2(xx[,col], na.rm=na.rm),
                      mean = mean   (xx[,col], na.rm=na.rm),
                      sd   = sd     (xx[,col], na.rm=na.rm)
                   )
                 },
                 measurevar,
                 na.rm
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean"=measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


allIHC <- read.csv(file.choose(),header=T)
distr <- subset(allIHC,Noise == '90dBNE' | Noise == '96dBSA' | Noise == '96dBZO'| Noise == 'NoNoise')

RibAbneuralPerSE <- summarySE(distr,measurevar='RibbonAbneuralPercent',groupvars=c('Noise','Frequency'))
ribAbneuralPerPlot <- ggplot(RibAbneuralPerSE, aes(x=Frequency, y=RibbonAbneuralPercent, color=Noise)) +
  geom_errorbar(aes(ymin=RibbonAbneuralPercent-se, ymax=RibbonAbneuralPercent+se), width=.3)+
  geom_line() +
  geom_point(size=3, shape=21, fill='white') +
  xlab('Frequency (kHz)') +
  ylab('Percentage of abneural side ribbon') +
  scale_colour_hue()+
  theme_bw() +
  theme(legend.justification=c(0,1),legend.position='top')

RibNeuralPerSE <- summarySE(distr,measurevar='RibbonNeuralPercent',groupvars=c('Noise','Frequency'))
ribNeuralPerPlot <- ggplot(RibNeuralPerSE, aes(x=Frequency, y=RibbonNeuralPercent, color=Noise)) +
  geom_errorbar(aes(ymin=RibbonNeuralPercent-se, ymax=RibbonNeuralPercent+se), width=.3)+
  geom_line() +
  geom_point(size=3, shape=21, fill='white') +
  xlab('Frequency (kHz)') +
  ylab('Percentage of neural side ribbon') +
  scale_colour_hue()+
  theme_bw() +
  theme(legend.justification=c(0,1),legend.position='top')

RibApicalPerSE <- summarySE(distr,measurevar='RibbonApicalPercent',groupvars=c('Noise','Frequency'))
ribApicalPerPlot <- ggplot(RibApicalPerSE, aes(x=Frequency, y=RibbonApicalPercent, color=Noise)) +
  geom_errorbar(aes(ymin=RibbonApicalPercent-se, ymax=RibbonApicalPercent+se), width=.3)+
  geom_line() +
  geom_point(size=3, shape=21, fill='white') +
  xlab('Frequency (kHz)') +
  ylab('Percentage of apical side ribbon') +
  scale_colour_hue()+
  theme_bw() +
  theme(legend.justification=c(0,1),legend.position='top')

RibBasalPerSE <- summarySE(distr,measurevar='RibbonBasalPercent',groupvars=c('Noise','Frequency'))
ribBasalPerPlot <- ggplot(RibBasalPerSE, aes(x=Frequency, y=RibbonBasalPercent, color=Noise)) +
  geom_errorbar(aes(ymin=RibbonBasalPercent-se, ymax=RibbonBasalPercent+se), width=.3)+
  geom_line() +
  geom_point(size=3, shape=21, fill='white') +
  xlab('Frequency (kHz)') +
  ylab('Percentage of Basal side ribbon') +
  scale_colour_hue()+
  theme_bw() +
  theme(legend.justification=c(0,1),legend.position='top')

multiplot(ribApicalPerPlot, ribAbneuralPerPlot,ribBasalPerPlot, ribNeuralPerPlot, cols=2)

RecNeural2PerSE <- summarySE(distr,measurevar='ReceptorNeural2Percent',groupvars=c('Noise','Frequency'))
recNeural2PerPlot <- ggplot(RecNeural2PerSE, aes(x=Frequency, y=ReceptorNeural2Percent, color=Noise)) +
  geom_errorbar(aes(ymin=ReceptorNeural2Percent-se, ymax=ReceptorNeural2Percent+se), width=.3)+
  geom_line() +
  geom_point(size=3, shape=21, fill='white') +
  xlab('Frequency (kHz)') +
  ylab('Percentage of Neural side receptor') +
  scale_colour_hue()+
  theme_bw() +
  theme(legend.justification=c(0,1),legend.position='top')

RibNeural2PerSE <- summarySE(distr,measurevar='RibbonNeural2Percent',groupvars=c('Noise','Frequency'))
ribNeural2PerPlot <- ggplot(RibNeural2PerSE, aes(x=Frequency, y=RibbonNeural2Percent, color=Noise)) +
  geom_errorbar(aes(ymin=RibbonNeural2Percent-se, ymax=RibbonNeural2Percent+se), width=.3)+
  geom_line() +
  geom_point(size=3, shape=21, fill='white') +
  xlab('Frequency (kHz)') +
  ylab('Percentage of neural side ribbon') +
  scale_colour_hue()+
  theme_bw() +
  theme(legend.justification=c(0,1),legend.position='top')

multiplot( ribNeural2PerPlot, recNeural2PerPlot,cols=2)

distr$Frequency <- as.factor(distr$Frequency)
theAOV <- aov(RibbonAbneuralPercent ~ Noise * Frequency, data=distr)
ribAbneuralPerANOVA <- anova(theAOV)
ribAbneuralPerTukey <- TukeyHSD(theAOV,'Noise')

theAOV <- aov(RibbonNeuralPercent ~ Noise * Frequency, data=distr)
ribNeuralPerANOVA <- anova(theAOV)
ribNeuralPerTukey <- TukeyHSD(theAOV,'Noise')

theAOV <- aov(RibbonApicalPercent ~ Noise * Frequency, data=distr)
ribApicalPerANOVA <- anova(theAOV)
ribApicalPerTukey <- TukeyHSD(theAOV,'Noise')

theAOV <- aov(RibbonBasalPercent ~ Noise * Frequency, data=distr)
ribBasalPerANOVA <- anova(theAOV)
ribBasalPerTukey <- TukeyHSD(theAOV,'Noise')
