library('ggplot2')
allIHC <- read.csv(file.choose(),header=T,sep='\t')
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
