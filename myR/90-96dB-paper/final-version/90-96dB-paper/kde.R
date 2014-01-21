plot(density(subset(the20k, Noise == 'NoNoise',select=c('MinDist'))$MinDist),col='red',xlab='Minimal Distance (um)', 
     ylab='Density', main=NULL, xlim=c(0,8))
lines(density(subset(the20k, Noise == '90dBNE',select=c('MinDist'))$MinDist),col='green',pch='.')
lines(density(subset(the20k, Noise == '96dBSA',select=c('MinDist'))$MinDist),col='steelblue',pch='.')
lines(density(subset(the20k, Noise == '96dBZO',select=c('MinDist'))$MinDist),col='purple',pch='.')
