# Automatic generate folder
theRnw <- 'LabJournal2014.Rnw'
monthpath <- format(Sys.time(), "%Y-%m")
daypath <- format(Sys.time(), "%Y-%m-%d-%a")
allmonth <- list.dirs(full.names=F, recursive=F)
if (! monthpath %in% allmonth) {
  system(paste('mkdir', monthpath))
}
daypath <- file.path(monthpath, daypath)
system(paste('mkdir', daypath))
system(paste('cp', 'Journal.tex', daypath))
jrnl <- file.path('.', daypath, 'Journal.tex')
nlines <- readLines(jrnl, warn=F)
nlines[1] <- paste("\\labday{", format(Sys.time(), "%A, %d %B %Y"), "}", sep='')
writeLines(nlines,jrnl)
newInsert <- substr(jrnl, 1, nchar(jrnl)-4)
nlines <- readLines(theRnw)
pos <- max(grep('include',nlines))
if (pos < 0) {
  pos <- grep('mainmatter',nlines)+1
}
if (grep(newInsert,nlines) < 0) {
  nlines <- c(nlines[1:pos], paste("\\include{",newInsert,"}", sep=''),
             nlines[(pos+1):length(nlines)])
  writeLines(nlines,theRnw)
}