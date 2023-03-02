path = '/mnt/ExtraDrive1/Work/desktop_data/2022_papers/hglmm/'
setwd(path)

library(classInt)
library(viridis)
library(sf)
library(xtable)
source(paste0(path,'autocorr_functions.R'))
source(paste0(path,'logLik_Laplace.R'))
source(paste0(path,'estpred.R'))
source(paste0(path,'addBreakColorLegend.R'))

load(paste0('seal_counts.rda'))

y = seal_counts$count
X = model.matrix(~ polyid + I(time_from_low/60) + I((time_from_low/60)^2) + 
	hrstd + I(hrstd^2), data = seal_counts)
#Z1 = model.matrix(~ -1 + polyid, data = DF)
Z1 = model.matrix(~ - 1 + polyid:yrfact, data = seal_counts)
mask = model.matrix(~ -1 + polyid, data = seal_counts)
mask = mask %*% t(mask)
distmat = as.matrix(dist(seal_counts[, 'yr']))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  Poisson Regression
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

maxvar = 2*var(log(y+1))
maxrange = 0.999

source(paste0(path,'autocorr_functions.R'))
source(paste0(path,'logLik_Laplace.R'))
theta = rep(-2, times = 4)
#undebug(logLik_Laplace)
optout_poisson = optim(theta, logLik_Laplace, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_AR1, maxvar = maxvar, maxrange = maxrange, 
	maxphi = maxphi, Z1 = Z1, mask = mask, 
	family = 'poisson', mlmeth = 'reml')
optout_poisson$value

#undebug(estpred)
epout_poisson = estpred(optout_poisson$par, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, 
	maxphi = maxphi, Z1 = Z1, mask = mask, 
	family = 'poisson')

epout_poisson$theta_trans

FEtable_poisson = cbind(epout_poisson$betahat, 
	sqrt(diag(epout_poisson$covbeta)), 
	sqrt(diag(epout_poisson$covbeta_adj)),
	abs(epout_poisson$betahat/sqrt(diag(epout_poisson$covbeta_adj))),
	2*(1-pt(abs(epout_poisson$betahat/sqrt(diag(epout_poisson$covbeta_adj))), 
		df = length(y) - length(X[1,]))))

FEtable_poisson[75:78,]

print(
	xtable(
		FEtable_poisson[75:78,],
     align = c('l',rep('l', times = 5)),
      digits = c(0,rep(3, times = 4),4),
      caption = 'Aymptotic correlation matrix for covariance parameters',
      label = 'tab:AsySE'
    ),
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

min(DF$hrstd)
x = (-15:15)/10
plot(x,  - 1.153*x - 1.498*x^2)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  Negative Binomial Regression
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

maxvar = 2*var(log(y+1))
maxrange = 0.999
maxphi = 100

source(paste0(path,'autocorr_functions.R'))
source(paste0(path,'logLik_Laplace.R'))
theta = rep(-2, times = 4)
#undebug(logLik_Laplace)
optout_negbinomial = optim(theta, logLik_Laplace, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_AR1, maxvar = maxvar, maxrange = maxrange, 
	maxphi = maxphi, Z1 = Z1, mask = mask, 
	family = 'negbinomial', mlmeth = 'reml', use.nugget = FALSE)
optout_negbinomial$value

#undebug(estpred)
epout_negbinomial = estpred(optout_negbinomial$par, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, 
	maxphi = maxphi, Z1 = Z1, mask = mask, 
	family = 'negbinomial', use.nugget = FALSE)

epout_negbinomial$theta_trans

FEtable_negbin = cbind(epout_negbinomial$betahat, 
	sqrt(diag(epout_negbinomial$covbeta)), 
	sqrt(diag(epout_negbinomial$covbeta_adj)),
	abs(epout_negbinomial$betahat/sqrt(diag(epout_negbinomial$covbeta_adj))),
	2*(1-pt(abs(epout_negbinomial$betahat/
		sqrt(diag(epout_negbinomial$covbeta_adj))), 
		df = length(y) - length(X[1,]))))
		
FEtable_negbin[75:78,]

print(
	xtable(
		FEtable_negbin[75:78,],
     align = c('l',rep('l', times = 5)),
      digits = c(0,3,5,3,3,4),
      caption = 'Aymptotic correlation matrix for covariance parameters',
      label = 'tab:AsySE'
    ),
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)
	
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  Plot Fitted Covariate Effects
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = "figures/seals_explanvar"
pdf(paste0(file_name,'.pdf'), width = 13, height = 7)

	layout(matrix(1:2, nrow = 1, byrow = TRUE))

	min(DF$hrstd)
	x_hrstd = (-150:150)/100
	par(mar = c(5,5,5,1))
	plot(x_hrstd,  epout_poisson$betahat[77]*x_hrstd + 
		epout_negbinomial$betahat[78]*x_hrstd^2, type = 'l', lwd = 3,
		xlab = 'Hour of Day (Standardized)',
		ylab = 'Log of Relative Probability',
		cex.lab = 2, cex.axis = 1.5)
	mtext('A', cex = 4, side = 3, adj = -.2, padj = -.3)
		
	x_timelow = (-200:200)/100
	par(mar = c(5,5,5,1))
	plot(x_timelow,  epout_poisson$betahat[75]*x_timelow + 
		epout_negbinomial$betahat[76]*x_timelow^2, type = 'l', lwd = 3,
		xlab = 'Time from Low Tide (Minutes)',
		ylab = 'Log of Relative Probability',
		cex.lab = 2, cex.axis = 1.5)
	mtext('B', cex = 4, side = 3, adj = -.2, padj = -.3)

	layout(1)

dev.off()

system(paste0('pdfcrop ','\'',
	path,file_name,'.pdf','\''))
system(paste0('cp ','\'',
	path,file_name,'-crop.pdf','\' ','\'',
	path,file_name,'.pdf','\''))
system(paste0('rm ','\'',
		path,file_name,'-crop.pdf','\''))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  Plot Fitted Time Series
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = "figures/seals_predw"
pdf(paste0(file_name,'.pdf'), width = 13, height = 13)

layout(matrix(1:4, nrow = 2, byrow = TRUE))

# plot 1

polyid_no = 8
DFp = seal_counts[1:19,]
DFp$yr = 1998:2016
DFp$time_from_low = 0
DFp$hrstd = 0
DFp$polyid = levels(seal_counts$polyid)[polyid_no]
DFp$count = NA

Xp = model.matrix(~ polyid + time_from_low + I(time_from_low^2) + hrstd + I(hrstd^2), 
	data = rbind(DFp,seal_counts))[1:19,]
Z1p = model.matrix(~ - 1 + polyid:yrfact, data = rbind(DFp,seal_counts))[1:19,]
dist_op = as.matrix(dist(matrix(rbind(seal_counts,DFp)$yr, ncol = 1)))[1:length(y),
	(length(y) + 1):(length(y) + 19)]
dist_pp = as.matrix(dist(matrix(DFp$yr, ncol = 1)))
maskp = outer((seal_counts$polyid == levels(seal_counts$polyid)[polyid_no])*1, 
	rep(1, times = 19))

#undebug(estpred)
epout_poisson = estpred(optout_poisson$par, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, 
	maxphi = maxphi, Z1 = Z1, mask = mask, 
	family = 'poisson', Xp = Xp, Z1p = Z1p, Z2p = Z2p, 
	dist_op = dist_op, dist_pp = dist_pp, maskp = maskp)
cbind(epout_poisson$w_pred, epout_poisson$wpred_se, 
	epout_poisson$wpred_se_adj)

pid_sub = seal_counts[seal_counts$polyid == 
	levels(seal_counts$polyid)[polyid_no], c('yr','count')]
UB = exp(epout_poisson$w_pred + 1.645*epout_poisson$wpred_se_adj)
LB = exp(epout_poisson$w_pred - 1.645*epout_poisson$wpred_se_adj)
par(mar = c(5,5,5,1))
plot(pid_sub,  xlim = c(1998, 2016), cex = 3, cex.lab = 2, cex.axis = 1.5,
	ylim = c(min(LB), max(UB)),
	main = levels(seal_counts$polyid)[polyid_no], cex.main = 2,
	xlab = 'Year')
points(1998:2016, exp(epout_poisson$w_pred), cex = 2, pch = 19)
lines(1998:2016, exp(epout_poisson$w_pred), lwd = 2)
lines(1998:2016, UB, lwd = 2, lty = 2)
lines(1998:2016, LB, lwd = 2, lty = 2)

# plot 2

polyid_no = 9
DFp = seal_counts[1:19,]
DFp$yr = 1998:2016
DFp$time_from_low = 0
DFp$hrstd = 0
DFp$polyid = levels(seal_counts$polyid)[polyid_no]
DFp$count = NA

Xp = model.matrix(~ polyid + time_from_low + I(time_from_low^2) + hrstd + I(hrstd^2), 
	data = rbind(DFp,seal_counts))[1:19,]
Z1p = model.matrix(~ - 1 + polyid:yrfact, data = rbind(DFp,seal_counts))[1:19,]
dist_op = as.matrix(dist(matrix(rbind(seal_counts,DFp)$yr, ncol = 1)))[1:length(y),
	(length(y) + 1):(length(y) + 19)]
dist_pp = as.matrix(dist(matrix(DFp$yr, ncol = 1)))
maskp = outer((seal_counts$polyid == levels(seal_counts$polyid)[polyid_no])*1, 
	rep(1, times = 19))

#undebug(estpred)
epout_poisson = estpred(optout_poisson$par, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, 
	maxphi = maxphi, Z1 = Z1, mask = mask, 
	family = 'poisson', Xp = Xp, Z1p = Z1p, Z2p = Z2p, 
	dist_op = dist_op, dist_pp = dist_pp, maskp = maskp)
cbind(epout_poisson$w_pred, epout_poisson$wpred_se, 
	epout_poisson$wpred_se_adj)

pid_sub = seal_counts[seal_counts$polyid == 
	levels(seal_counts$polyid)[polyid_no], c('yr','count')]
UB = exp(epout_poisson$w_pred + 1.645*epout_poisson$wpred_se_adj)
LB = exp(epout_poisson$w_pred - 1.645*epout_poisson$wpred_se_adj)
par(mar = c(5,5,5,1))
plot(pid_sub,  xlim = c(1998, 2016), cex = 3, cex.lab = 2, cex.axis = 1.5,
	ylim = c(min(LB), max(UB)),
	main = levels(seal_counts$polyid)[polyid_no], cex.main = 2,
	xlab = 'Year')
points(1998:2016, exp(epout_poisson$w_pred), cex = 2, pch = 19)
lines(1998:2016, exp(epout_poisson$w_pred), lwd = 2)
lines(1998:2016, UB, lwd = 2, lty = 2)
lines(1998:2016, LB, lwd = 2, lty = 2)

# plot 3

polyid_no = 15
DFp = seal_counts[1:19,]
DFp$yr = 1998:2016
DFp$time_from_low = 0
DFp$hrstd = 0
DFp$polyid = levels(seal_counts$polyid)[polyid_no]
DFp$count = NA

Xp = model.matrix(~ polyid + time_from_low + I(time_from_low^2) + hrstd + I(hrstd^2), 
	data = rbind(DFp,seal_counts))[1:19,]
Z1p = model.matrix(~ - 1 + polyid:yrfact, data = rbind(DFp,seal_counts))[1:19,]
dist_op = as.matrix(dist(matrix(rbind(seal_counts,DFp)$yr, ncol = 1)))[1:length(y),
	(length(y) + 1):(length(y) + 19)]
dist_pp = as.matrix(dist(matrix(DFp$yr, ncol = 1)))
maskp = outer((seal_counts$polyid == levels(seal_counts$polyid)[polyid_no])*1, 
	rep(1, times = 19))

#undebug(estpred)
epout_poisson = estpred(optout_poisson$par, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, 
	maxphi = maxphi, Z1 = Z1, mask = mask, 
	family = 'poisson', Xp = Xp, Z1p = Z1p, Z2p = Z2p, 
	dist_op = dist_op, dist_pp = dist_pp, maskp = maskp)
cbind(epout_poisson$w_pred, epout_poisson$wpred_se, 
	epout_poisson$wpred_se_adj)

pid_sub = seal_counts[seal_counts$polyid == 
	levels(seal_counts$polyid)[polyid_no], c('yr','count')]
UB = exp(epout_poisson$w_pred + 1.645*epout_poisson$wpred_se_adj)
LB = exp(epout_poisson$w_pred - 1.645*epout_poisson$wpred_se_adj)
par(mar = c(5,5,5,1))
plot(pid_sub,  xlim = c(1998, 2016), cex = 3, cex.lab = 2, cex.axis = 1.5,
	ylim = c(min(LB), max(UB)),
	main = levels(seal_counts$polyid)[polyid_no], cex.main = 2,
	xlab = 'Year')
points(1998:2016, exp(epout_poisson$w_pred), cex = 2, pch = 19)
lines(1998:2016, exp(epout_poisson$w_pred), lwd = 2)
lines(1998:2016, UB, lwd = 2, lty = 2)
lines(1998:2016, LB, lwd = 2, lty = 2)

# plot 4

polyid_no = 16
DFp = seal_counts[1:19,]
DFp$yr = 1998:2016
DFp$time_from_low = 0
DFp$hrstd = 0
DFp$polyid = levels(seal_counts$polyid)[polyid_no]
DFp$count = NA

Xp = model.matrix(~ polyid + time_from_low + I(time_from_low^2) + hrstd + I(hrstd^2), 
	data = rbind(DFp,seal_counts))[1:19,]
Z1p = model.matrix(~ - 1 + polyid:yrfact, data = rbind(DFp,seal_counts))[1:19,]
dist_op = as.matrix(dist(matrix(rbind(seal_counts,DFp)$yr, ncol = 1)))[1:length(y),
	(length(y) + 1):(length(y) + 19)]
dist_pp = as.matrix(dist(matrix(DFp$yr, ncol = 1)))
maskp = outer((seal_counts$polyid == levels(seal_counts$polyid)[polyid_no])*1, 
	rep(1, times = 19))

#undebug(estpred)
epout_poisson = estpred(optout_poisson$par, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, 
	maxphi = maxphi, Z1 = Z1, mask = mask, 
	family = 'poisson', Xp = Xp, Z1p = Z1p, Z2p = Z2p, 
	dist_op = dist_op, dist_pp = dist_pp, maskp = maskp)
cbind(epout_poisson$w_pred, epout_poisson$wpred_se, 
	epout_poisson$wpred_se_adj)

pid_sub = seal_counts[seal_counts$polyid == 
	levels(seal_counts$polyid)[polyid_no], c('yr','count')]
UB = exp(epout_poisson$w_pred + 1.645*epout_poisson$wpred_se_adj)
LB = exp(epout_poisson$w_pred - 1.645*epout_poisson$wpred_se_adj)
par(mar = c(5,5,5,1))
plot(pid_sub,  xlim = c(1998, 2016), cex = 3, cex.lab = 2, cex.axis = 1.5,
	ylim = c(min(LB), max(UB)),
	main = levels(seal_counts$polyid)[polyid_no], cex.main = 2,
	xlab = 'Year')
points(1998:2016, exp(epout_poisson$w_pred), cex = 2, pch = 19)
lines(1998:2016, exp(epout_poisson$w_pred), lwd = 2)
lines(1998:2016, UB, lwd = 2, lty = 2)
lines(1998:2016, LB, lwd = 2, lty = 2)

	layout(1)

dev.off()

system(paste0('pdfcrop ','\'',
	path,file_name,'.pdf','\''))
system(paste0('cp ','\'',
	path,file_name,'-crop.pdf','\' ','\'',
	path,file_name,'.pdf','\''))
system(paste0('rm ','\'',
		path,file_name,'-crop.pdf','\''))
