library(HSsurv2018)
library(unitted)
library(lubridate)

basedir1 = '/media/jay/ExtraDrive1/00NMML/activePapers/HSsurv2018'
basedir2 = '/HSsurv2018_package/HSsurv2018'
basedir = paste0(basedir1, basedir2)
data(dHOterr)
table(dHOterr$stockid)
data(dterr)

table(dterr$stockid)
dstk = dterr[dterr$stockid == 11,]
unique(as.character(dstk$polyid))
table(as.character(dstk$polyid))
#ind_one = which(table(as.character(dstk$polyid)) == 1)
#dstk = dstk[!dstk$polyid %in% names(ind_one),]

table(as.character(dstk$yr))
#dstk[dstk$polyid == 'JF00','yr']
subset = c('JE01', 'JE02', 'JE03','JE04','JE05','JE06','JE07','JE08',
	'JE20', 'JE21', 'JE22','JE23','JF00','JF01','JF02','JF03','JF04',
	'JF05','JF06','JF07','JF08','JF09','JF10','JF11','JF12','JF13')
#DF = dstk[dstk$polyid %in% subset,c('polyid','count',
#	'time_from_low','yr','day','hr','daystd','hrstd')]
DF = dstk[,c('polyid','count',
	'time_from_low','yr','day','hr','daystd','hrstd')]
DF$polyid = as.factor(as.character(DF$polyid))
DF$yrfact = as.factor(as.character(DF$yr))

seal_counts = DF

save(seal_counts, file = 'seal_counts.rda')
