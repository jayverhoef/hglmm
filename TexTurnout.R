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

load(paste0('Tex_turnout80.rda'))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  Scatterplot of variables
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = "figures/TexTurn_scatter"
pdf(paste0(file_name,'.pdf'), width = 12, height = 4.5)

cex_all = 1.5
cexlab_all = 2.2
cexaxis_all = 1.7
layout(matrix(1:3, nrow = 1), widths = c(1.2,1,1))
	par(mar = c(5,5,0,0))
	plot(Tex_turnout80$pc_college, logit(Tex_turnout80$pc_turnout), pch = 19,
		ylab = 'Logit(Proportion Turnout)',
		xlab = 'Proportion with College Degree',
		cex.lab = cexlab_all, cex.axis = cexaxis_all, cex = cex_all)
	par(mar = c(5,0,0,0))
	plot(Tex_turnout80$pc_homeownership^3, logit(Tex_turnout80$pc_turnout),
		pch = 19, ylab = '', yaxt = 'n',
		xlab = '(Proportion Home Ownership)^3',
		cex.lab = cexlab_all, cex.axis = cexaxis_all, cex = cex_all)
	plot(log(Tex_turnout80$pc_income), logit(Tex_turnout80$pc_turnout),
		pch = 19, ylab = '', yaxt = 'n',
		xlab = 'log(Per Capita Income)',
		cex.lab = cexlab_all, cex.axis = cexaxis_all, cex = cex_all)
layout(1)

dev.off()

system(paste0('pdfcrop ','\'',
	path,file_name,'.pdf','\''))
system(paste0('cp ','\'',
	path,file_name,'-crop.pdf','\' ','\'',
	path,file_name,'.pdf','\''))
system(paste0('rm ','\'',
		path,file_name,'-crop.pdf','\''))

xy = st_coordinates(Tex_turnout80)/1000


Distmat = as.matrix(dist(xy))
# create first-order neighbor matrix (rook's move) from distances
W = (Distmat < 150)*1
diag(W) = 0
min(apply(W, 1, sum))
max(apply(W, 1, sum))

# matrices for CAR and SAR models, using row-standardization 
Mbar = diag(1/apply(W, 1, sum))
Wbar = W/apply(W,1,sum)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  CAR model for Bernoulli data
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

y = Tex_turnout80$bin_turnout
X = model.matrix(~ pc_college + I(pc_homeownership^3) + 
	I(log(Tex_turnout80$pc_income)), data = Tex_turnout80)

maxvar = 4*var(logit(Tex_turnout80$pc_turnout))
maxrange = 0.999
maxphi = 10000

theta = rep(0, times = 2)
source(paste0(path,'autocorr_functions.R'))
source(paste0(path,'logLik_Laplace.R'))
#undebug(logLik_Laplace)
optout_binCAR = optim(theta, logLik_Laplace, y = y, X = X, 
	autocor_fun = rho_CAR, maxvar = maxvar, maxrange = maxrange,
	family = 'binomial', mlmeth = 'reml', W = Wbar, M = Mbar,
	use.nugget = FALSE)
optout_binCAR$value


source(paste0(path,'estpred.R'))
#undebug(estpred)
epout_binCAR = estpred(optout_binCAR$par, y = y, X = X, 
	autocor_fun = rho_CAR, maxvar = maxvar, maxrange = maxrange, 
	family = 'binomial', W = Wbar, M = Mbar, use.nugget = FALSE)

epout_binCAR$theta_trans

FEtable_binCAR = cbind(epout_binCAR$betahat, 
	sqrt(diag(epout_binCAR$covbeta)), 
	sqrt(diag(epout_binCAR$covbeta_adj)),
	abs(epout_binCAR$betahat/sqrt(diag(epout_binCAR$covbeta_adj))),
	2*(1-pt(abs(epout_binCAR$betahat/sqrt(diag(epout_binCAR$covbeta_adj))), 
		df = length(y) - length(X[1,]))))
	

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  SAR model for Bernoulli data
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

y = Tex_turnout80$bin_turnout
X = model.matrix(~ pc_college + I(pc_homeownership^3) + 
	I(log(Tex_turnout80$pc_income)), data = Tex_turnout80)

maxvar = 4*var(logit(Tex_turnout80$pc_turnout))
maxrange = 0.999
maxphi = 10000

theta = rep(0, times = 2)
source(paste0(path,'autocorr_functions.R'))
source(paste0(path,'logLik_Laplace.R'))
#undebug(logLik_Laplace)
optout_binSAR = optim(theta, logLik_Laplace, y = y, X = X, 
	autocor_fun = rho_SAR, maxvar = maxvar, maxrange = maxrange,
	family = 'binomial', mlmeth = 'reml', W = Wbar, M = Mbar,
	use.nugget = FALSE)
optout_binSAR$value


source(paste0(path,'estpred.R'))
#undebug(estpred)
epout_binSAR = estpred(optout_binSAR$par, y = y, X = X, 
	autocor_fun = rho_SAR, maxvar = maxvar, maxrange = maxrange, 
	family = 'binomial', W = Wbar, M = Mbar, use.nugget = FALSE)

epout_binSAR$theta_trans

FEtable_binSAR = cbind(epout_binSAR$betahat, 
	sqrt(diag(epout_binSAR$covbeta)), 
	sqrt(diag(epout_binSAR$covbeta_adj)),
	abs(epout_binSAR$betahat/sqrt(diag(epout_binSAR$covbeta_adj))),
	2*(1-pt(abs(epout_binSAR$betahat/sqrt(diag(epout_binSAR$covbeta_adj))), 
		df = length(y) - length(X[1,]))))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  SAR model for proportion data
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

y = Tex_turnout80$pc_turnout
X = model.matrix(~ 1, data = Tex_turnout80)
X = model.matrix(~ pc_college + I(pc_homeownership^3) + 
	I(log(Tex_turnout80$pc_income)), data = Tex_turnout80)

maxvar = 2*var(logit(Tex_turnout80$pc_turnout))
maxrange = 0.999
maxphi = 10000

theta = rep(-2, times = 3)
source(paste0(path,'autocorr_functions.R'))
source(paste0(path,'logLik_Laplace.R'))
#undebug(logLik_Laplace)
optout_turnoutSAR = optim(theta, logLik_Laplace, y = y, X = X, 
	autocor_fun = rho_SAR, maxvar = maxvar, maxrange = maxrange, maxphi = maxphi,
	family = 'beta', mlmeth = 'reml', W = Wbar, M = Mbar,
	use.nugget = FALSE)
optout_turnoutSAR$value

source(paste0(path,'estpred.R'))
#undebug(estpred)
epout_turnoutSAR = estpred(optout_turnoutSAR$par, y = y, X = X, 
	autocor_fun = rho_SAR, maxvar = maxvar, maxrange = maxrange, maxphi = maxphi,
	family = 'beta', W = Wbar, M = Mbar, use.nugget = FALSE)
	
epout_turnoutSAR$theta_trans

FEtable_turnoutSAR = cbind(epout_turnoutSAR$betahat, 
	sqrt(diag(epout_turnoutSAR$covbeta)), 
	sqrt(diag(epout_turnoutSAR$covbeta_adj)),
	abs(epout_turnoutSAR$betahat/sqrt(diag(epout_turnoutSAR$covbeta_adj))),
	2*(1-pt(abs(epout_turnoutSAR$betahat/sqrt(diag(epout_turnoutSAR$covbeta_adj))), 
		df = length(y) - length(X[1,]))))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  CAR model for proportion data
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

y = Tex_turnout80$pc_turnout
X = model.matrix(~ 1, data = Tex_turnout80)
X = model.matrix(~ pc_college + I(pc_homeownership^3) + 
	I(log(Tex_turnout80$pc_income)), data = Tex_turnout80)

maxvar = 2*var(logit(Tex_turnout80$pc_turnout))
maxrange = 0.999
maxphi = 10000

theta = rep(0, times = 3)
source(paste0(path,'autocorr_functions.R'))
source(paste0(path,'logLik_Laplace.R'))
#undebug(logLik_Laplace)
optout_turnoutCAR = optim(theta, logLik_Laplace, y = y, X = X, 
	autocor_fun = rho_CAR, maxvar = maxvar, maxrange = maxrange, maxphi = maxphi,
	family = 'beta', mlmeth = 'reml', W = Wbar, M = Mbar,
	use.nugget = FALSE)
optout_turnoutCAR$value

source(paste0(path,'estpred.R'))
#undebug(estpred)
epout_turnoutCAR = estpred(optout_turnoutCAR$par, y = y, X = X, 
	autocor_fun = rho_CAR, maxvar = maxvar, maxrange = maxrange, maxphi = maxphi,
	family = 'beta', W = Wbar, M = Mbar, use.nugget = FALSE)
	
epout_turnoutCAR$theta_trans

FEtable_turnoutCAR = cbind(epout_turnoutCAR$betahat, 
	sqrt(diag(epout_turnoutCAR$covbeta)), 
	sqrt(diag(epout_turnoutCAR$covbeta_adj)),
	abs(epout_turnoutCAR$betahat/sqrt(diag(epout_turnoutCAR$covbeta_adj))),
	2*(1-pt(abs(epout_turnoutCAR$betahat/sqrt(diag(epout_turnoutCAR$covbeta_adj))), 
		df = length(y) - length(X[1,]))))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  Print Fixed Effects Tables
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

print(
	xtable(
		FEtable_binCAR,
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

print(
	xtable(
		FEtable_binSAR,
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

print(
	xtable(
		FEtable_turnoutCAR,
     align = c('l',rep('l', times = 5)),
      digits = c(0,rep(3, times = 4),4)
    ),
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

print(
	xtable(
		FEtable_turnoutSAR,
     align = c('l',rep('l', times = 5)),
      digits = c(0,rep(3, times = 4),4)
    ),
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  Maps of W-values
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = "figures/TexTurn_maps"
pdf(paste0(file_name,'.pdf'), width = 8.5, height = 11)

layout(matrix(1:12, ncol = 4, byrow = TRUE), widths = c(4,1,4,1))

brks_cex = 1.5
pts_cex = 2
leg_right = 0.5
lab_cex = 4
mar1 = c(0,2,0,2)
mar2 = c(0,0,0,0)
ybot = .2
ytop = .8

# A

par(mar = mar1)
plot(st_geometry(Tex_turnout80), type = 'n')
plot(st_geometry(Tex_turnout80[Tex_turnout80$bin == 0,]), add = TRUE,
	cex = pts_cex)
plot(st_geometry(Tex_turnout80[Tex_turnout80$bin == 1,]), add = TRUE,
	pch = 19, cex = pts_cex)
text(-961569, 1503770, 'A', cex = lab_cex, adj = c(0,1))
par(mar = mar2)
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
#addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
#  breaks = c(0,1,2), colors = c('white','black'), cex = brks_cex, 
#  printFormat = "1.0")

# D

f10 = classIntervals(logit(Tex_turnout80$pc_turnout), n = 10, style = 'fisher')
pal = viridis(length(f10$brks) - 1)
f10Colours = findColours(f10, pal)
par(mar = mar1)
plot(st_geometry(Tex_turnout80), col = f10Colours, pch = 19, cex = 2)
text(-961569, 1503770, 'D', cex = lab_cex, adj = c(0,1))
par(mar = mar2)
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = ybot, xright = .2, ytop = ytop,
  breaks = f10$brks, colors = pal, cex = brks_cex, 
  printFormat = "1.2")

# B

par(mar = c(0,0,0,0))
f10 = classIntervals(epout_binSAR$w, n = 10, style = 'fisher')
pal = viridis(length(f10$brks) - 1)
f10Colours = findColours(f10, pal)
par(mar = mar1)
plot(st_geometry(Tex_turnout80), col = f10Colours, pch = 19, cex = pts_cex)
text(-961569, 1503770, 'B', cex = lab_cex, adj = c(0,1))
par(mar = mar2)
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = ybot, xright = .2, ytop = ytop,
  breaks = f10$brks, colors = pal, cex = brks_cex, 
  printFormat = "1.2")

# E
par(mar = c(0,0,0,0))
f10 = classIntervals(epout_turnoutSAR$w, n = 10, style = 'fisher')
pal = viridis(length(f10$brks) - 1)
f10Colours = findColours(f10, pal)
par(mar = mar1)
plot(st_geometry(Tex_turnout80), col = f10Colours, pch = 19, cex = pts_cex)
text(-961569, 1503770, 'E', cex = lab_cex, adj = c(0,1))
par(mar = mar2)
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = ybot, xright = .2, ytop = ytop,
  breaks = f10$brks, colors = pal, cex = brks_cex, 
  printFormat = "1.2")

#C

f10 = classIntervals(epout_binCAR$w, n = 10, style = 'fisher')
pal = viridis(length(f10$brks) - 1)
f10Colours = findColours(f10, pal)
par(mar = mar1)
plot(st_geometry(Tex_turnout80), col = f10Colours, pch = 19, cex = pts_cex)
text(-961569, 1503770, 'C', cex = lab_cex, adj = c(0,1))
par(mar = c(0,0,0,0))
par(mar = mar2)
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = ybot, xright = .2, ytop = ytop,
  breaks = f10$brks, colors = pal, cex = brks_cex, 
  printFormat = "1.2")


# F
f10 = classIntervals(epout_turnoutCAR$w, n = 10, style = 'fisher')
pal = viridis(length(f10$brks) - 1)
f10Colours = findColours(f10, pal)
par(mar = mar1)
plot(st_geometry(Tex_turnout80), col = f10Colours, pch = 19, cex = pts_cex)
text(-961569, 1503770, 'F', cex = lab_cex, adj = c(0,1))
par(mar = mar2)
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = ybot, xright = .2, ytop = ytop,
  breaks = f10$brks, colors = pal, cex = brks_cex, 
  printFormat = "1.2")

layout(1)

dev.off()

system(paste0('pdfcrop ','\'',
	path,file_name,'.pdf','\''))
system(paste0('cp ','\'',
	path,file_name,'-crop.pdf','\' ','\'',
	path,file_name,'.pdf','\''))
system(paste0('rm ','\'',
		path,file_name,'-crop.pdf','\''))


cbind(epout_binSAR$w, epout_binCAR$w, epout_turnoutSAR$w, epout_binCAR$w

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  Fitted beta distribution
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

store_mat = matrix(NA, ncol = length(y), nrow = 999)
phi = 48.7
for(i in 1:length(y)) store_mat[,i] = dbeta((1:999)/1000, 
	shape1 = phi*expit(epout_turnoutSAR$w[i]), 
	shape2 = phi*(1- expit(epout_turnoutSAR$w[i])))

file_name = "figures/TexTurn_histpdf"
pdf(paste0(file_name,'.pdf'), width = 6, height = 6)

par(mar = c(5,5,1,1))
hist(Tex_turnout80$pc_turnout, freq = FALSE, xlim = c(0,1),
	ylab = expression(paste("[y|", mu, ", ", phi,"] = Beta(",mu,", ",phi,")")), xlab = 'y', cex.lab = 2, cex.axis = 1.5, breaks = (0:10/10), ylim = c(0,8),
	main = '')
#lines((1:999)/1000, apply(store_mat,1,mean), lwd = 3)
lines((1:999)/1000, dbeta((1:999)/1000, shape1 = phi*expit(logit(.5)), 
	shape2 = phi*(1- expit(logit(.5)))), lwd = 3, lty = 1)
lines((1:999)/1000, dbeta((1:999)/1000, shape1 = phi*expit(logit(.2)), 
	shape2 = phi*(1- expit(logit(.3)))), lwd = 3, lty = 2)
lines((1:999)/1000, dbeta((1:999)/1000, shape1 = phi*expit(logit(.9)), 
	shape2 = phi*(1- expit(logit(.8)))), lwd = 3, lty = 3)

dev.off()

system(paste0('pdfcrop ','\'',
	path,file_name,'.pdf','\''))
system(paste0('cp ','\'',
	path,file_name,'-crop.pdf','\' ','\'',
	path,file_name,'.pdf','\''))
system(paste0('rm ','\'',
		path,file_name,'-crop.pdf','\''))

lines((1:999)/1000, dbeta((1:999)/1000, shape1 = phi*expit(0), 
	shape2 = phi*(1- expit(0))), lwd = 3, lty = 1)
lines((1:999)/1000, dbeta((1:999)/1000, shape1 = phi*expit(logit(.2)), 
	shape2 = phi*(1- expit(logit(.3)))), lwd = 3, lty = 2)
lines((1:999)/1000, dbeta((1:999)/1000, shape1 = phi*expit(logit(.9)), 
	shape2 = phi*(1- expit(logit(.8)))), lwd = 3, lty = 3)


