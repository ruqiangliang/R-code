library(tools)
setwd("/media/ruqiang/OS/Users/Ruqiang/Desktop/Toshiba-USB/LMD6K/2013-02-06")
cochlea <- dir('.')
for (i in 1: length(cochlea)){
  setwd(file.path(getwd(), cochlea[i], fsep= .Platform$file.sep))
  pic <- dir('.')
  for (ff in 1:length(pic)){
    theJpeg <- paste(file_path_sans_ext(pic[ff]),'.jpg',sep='')
    system(paste("convert", pic[ff], theJpeg ,sep=' '))
  }
  setwd('..')
}