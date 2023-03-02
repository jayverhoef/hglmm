path = '/mnt/ExtraDrive1/Work/desktop_data/2022_papers/hglmm/'
setwd(path)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                         Get the Data
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# attach data library
library(ZVHdata)
library(sf)
library(viridis)
library(classInt)
library(colorspace)
library(gstat)
library(spmodel)
library(stringr)
library(xtable)

# load data for graphics and analysis
data(MOSSobs)
data(MOSSpreds)
data(CAKRboundary)

# transform some of the variables
DF = data.frame(st_drop_geometry(MOSSobs), 
	easting = st_coordinates(MOSSobs)[,1]/1e+3,
	northing = st_coordinates(MOSSobs)[,2]/1e+3)
DF$year = as.factor(DF$year)
DF$field_dup = as.factor(DF$field_dup)
DF$dist2road = log(DF$dist2road)

sum(DF$year == '2001')
sum(DF$year == '2006')

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  Sample Design with Autocorrelation Points
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

DF01 = DF[DF$year == '2001',]
indAC = str_detect(DF01$sample, 'AC')

#file_name = 'figures/Moss_autocorr_samples'
#pdf(paste0(file_name,'.pdf'), width = 16, height = 8)

layout(matrix(1:2, nrow = 1))

	par(mar = c(0,0,1,0))
	plot(st_geometry(CAKRboundary))
	rect(-420000, 1985000, -410000, 1995000, lwd = 2, col = 'grey80')
	plot(st_geometry(MOSSobs[MOSSobs$year == 2001,]), add = TRUE, pch = 19)
	text(-437000, 2011000, 'A', cex = 4)

	par(mar = c(5,5,3,1))
	plot(DF01[!indAC, c('easting','northing')], pch = 3, cex = 3,
		xlim = c(-420,-410), ylim = c(1985,1995), cex.lab = 2, cex.axis = 1.5,
		xlab = 'Northing', ylab = 'Easting', lwd = 2)
	points(DF01[indAC,c('easting','northing')], pch = 4, cex = 3, lwd = 2)
	mtext('B', adj = -.2, cex = 4, padj = .4)

	layout(1)

#dev.off()

#system(paste0('pdfcrop ','\'',SLEDbook_path,
#	sec_path,file_name,'.pdf','\''))
#system(paste0('cp ','\'',SLEDbook_path,
#	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
#	sec_path,file_name,'.pdf','\''))
#system(paste0('rm ','\'',SLEDbook_path,
#		sec_path,file_name,'-crop.pdf','\''))

logit = function(p){log(p/(1-p))}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  log-transformation
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

PbOut = splm(I(log(Pb)) ~ year + dist2road + dist2road:sideroad, data = DF,
	xcoord = 'easting', ycoord = 'northing', spcov_type = 'exponential',
	partition_factor = ~ year,
	random = ~ (1 | sample) + (1 |sample:field_dup),
	control = list(reltol = 1e-12), estmethod = 'ml')
summary(PbOut)

fitted(PbOut)
vcov(PbOut)

X = model.matrix(~ year + dist2road + dist2road:sideroad, data = DF)
diag(X %*% vcov(PbOut) %*% t(X))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  Gamma Regression
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

y = DF$Pb
X = model.matrix(~ year + sideroad + dist2road + dist2road:sideroad, data = DF)
Z1 = model.matrix(~ -1 + sample, data = DF)
Z2 = model.matrix(~ - 1 + sample:field_dup, data = DF)
mask = model.matrix(~ -1 + I(as.factor(year)), data = DF)
mask = mask %*% t(mask)
distmat = as.matrix(dist(DF[, c('easting', 'northing')]))

maxvar = 4*var(log(y+1))
maxrange = 4*max(distmat)
maxphi = 10000

# undebug(logLik_Laplace)
source(paste0(path,'autocorr_functions.R'))
source(paste0(path,'logLik_Laplace.R'))
theta = rep(-4, times = 6)
optout_gamma = optim(theta, logLik_Laplace, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, 
	maxphi = maxphi, Z1 = Z1, Z2 = Z2, mask = mask, 
	family = 'gamma', mlmeth = 'ml')
optout_gamma$value

y = DF$Pb
X = model.matrix(~ year + dist2road + dist2road:sideroad, data = DF)
Z1 = model.matrix(~ -1 + sample, data = DF)
Z2 = model.matrix(~ - 1 + sample:field_dup, data = DF)
mask = model.matrix(~ -1 + I(as.factor(year)), data = DF)
mask = mask %*% t(mask)
distmat = as.matrix(dist(DF[, c('easting', 'northing')]))

maxvar = 4*var(log(y+1))
maxrange = 4*max(distmat)
maxphi = 10000

# undebug(logLik_Laplace)
source(paste0(path,'autocorr_functions.R'))
source(paste0(path,'logLik_Laplace.R'))
theta = rep(-4, times = 6)
optout_gamma = optim(theta, logLik_Laplace, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, 
	maxphi = maxphi, Z1 = Z1, Z2 = Z2, mask = mask, 
	family = 'gamma', mlmeth = 'ml')
optout_gamma$value

source(paste0(path,'estpred.R'))
epout_gamma = estpred(optout_gamma$par, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, 
	maxphi = maxphi, Z1 = Z1, Z2 = Z2, mask = mask, 
	family = 'gamma')

epout_gamma$theta_trans

cbind(epout_gamma$betahat, sqrt(diag(epout_gamma$covbeta)), 
sqrt(diag(epout_gamma$covbeta_adj)))

FEtable_gamma = cbind(epout_gamma$betahat, 
	sqrt(diag(epout_gamma$covbeta)), 
	sqrt(diag(epout_gamma$covbeta_adj)),
	abs(epout_gamma$betahat/sqrt(diag(epout_gamma$covbeta_adj))),
	2*(1-pt(abs(epout_gamma$betahat/
		sqrt(diag(epout_gamma$covbeta_adj))), 
		df = length(y) - length(X[1,]))))

print(
	xtable(
		FEtable_gamma,
     align = c('l',rep('l', times = 5)),
      digits = c(0,3,5,3,4,4),
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
#                  Inverse Gaussian Regression
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

maxvar = 4*var(log(y+1))
maxrange = 4*max(distmat)
maxphi = 10000

# undebug(logLik_Laplace)
source(paste0(path,'autocorr_functions.R'))
source(paste0(path,'logLik_Laplace.R'))
theta = rep(-4, times = 6)
optout_invgau = optim(theta, logLik_Laplace, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, 
	maxphi = maxphi, Z1 = Z1, Z2 = Z2, mask = mask, 
	family = 'invgauss', mlmeth = 'ml')
optout_invgau$value

source(paste0(path,'estpred.R'))
#optout$par[4] = 1
#undebug(estpred)
epout_invgau = estpred(optout_invgau$par, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, 
	maxphi = maxphi, Z1 = Z1, Z2 = Z2, mask = mask, 
	family = 'invgauss')

w_invgauss = epout_invgau$w

epout_invgau$theta_trans

FEtable_invgau = cbind(epout_invgau$betahat, 
	sqrt(diag(epout_invgau$covbeta)), 
	sqrt(diag(epout_invgau$covbeta_adj)),
	abs(epout_invgau$betahat/sqrt(diag(epout_invgau$covbeta_adj))),
	2*(1-pt(abs(epout_invgau$betahat/
		sqrt(diag(epout_invgau$covbeta_adj))), 
		df = length(y) - length(X[1,]))))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#    Distribution Comparisons
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = "figures/moss_densities"
pdf(paste0(file_name,'.pdf'), width = 13, height = 7)

	cex_lab = 2.5
	cex_axis = 1.8
	par_mar = c(5,5,5,1)
	layout(matrix(c(1,1,1,2,3,4), nrow = 2, byrow = TRUE))
	par(mar = par_mar)
	hist(y, xlim = c(0,1000), breaks = (0:25)*50, freq = FALSE, main = '',
		cex.lab = cex_lab, cex.axis = cex_axis, xlab = 'Concentration')
	mtext('A', cex = 4, side = 3, adj = -0.05, padj = -.2)
	yrang = (350:450)/10
	dinvgauss(y, mean = exp(w), shape = phi*exp(w), log = TRUE)
	plot(yrang, dgamma(yrang, shape = epout_gamma$theta_trans['phi'], 
		scale = 40/epout_gamma$theta_trans['phi']), type = 'l', lwd = 3,
		ylab = 'Density', cex.lab = cex_lab, cex.axis = cex_axis, 
		xlab = 'Concentration', ylim = c(0,.6), lty = 1)
	mtext('B', cex = 4, side = 3, adj = -0.19, padj = -.2)
	yrang = (1250:1750)/10
	plot(yrang, dgamma(yrang, shape = epout_gamma$theta_trans['phi'], 
		scale = 150/epout_gamma$theta_trans['phi']), type = 'l', lwd = 3,
		ylab = 'Density', cex.lab = cex_lab, cex.axis = cex_axis, 
		xlab = 'Concentration', ylim = c(0,.16), lty = 1)
	mtext('C', cex = 4, side = 3, adj = -0.19, padj = -.2)
	yrang = (7500:9500)/10
	plot(yrang, dgamma(yrang, shape = epout_gamma$theta_trans['phi'],
		scale = 850/epout_gamma$theta_trans['phi']), type = 'l', lty = 1,  lwd = 3,
		ylab = 'Density', cex.lab = cex_lab, cex.axis = cex_axis, 
		xlab = 'Concentration', ylim = c(0,0.031))
	mtext('D', cex = 4, side = 3, adj = -0.19, padj = -.2)

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
#    Predictions
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# prepare fixed quantities for predictions in strata
DFp = MOSSpreds
DFp$easting = st_coordinates(DFp)[,1]/1e+3
DFp$northing = st_coordinates(DFp)[,2]/1e+3
DFp = st_drop_geometry(DFp)
DFp$dist2road = log(DFp$dist2road)
DFp01 = DFp
DFp06 = DFp
DFp01[,'year'] = 2001
DFp06[,'year'] = 2006
DFpy = rbind(DFp01,DFp06)
DFpy$year = as.factor(DFpy$year)
	
#strata1
DFp1 = DFpy[DFpy$strat == 1,]
DF1 = rbind(DF[,c('easting','northing','year')],
	DFp1[,c('easting','northing','year')])
D1 = as.matrix(dist(DF1[,c('easting','northing')]))
part_mask1 = outer(DF1$year == '2001', DF1$year == '2001') + 
	outer(DF1$year == '2006', DF1$year == '2006')
n1 = dim(D1)[1]
distop1 = D1[1:length(y),(length(y) + 1):n1]
distpp1 = D1[(length(y) + 1):n1,(length(y) + 1):n1]
maskop1 = part_mask1[1:length(y),(length(y) + 1):n1]
maskpp1 = part_mask1[(length(y) + 1):n1,(length(y) + 1):n1]
Xp1 = model.matrix(~ year + dist2road + sideroad:dist2road, data = DFp1)

#strata2
DFp2 = DFpy[DFpy$strat == 2,]
DF2 = rbind(DF[,c('easting','northing','year')],
	DFp2[,c('easting','northing','year')])
D2 = as.matrix(dist(DF2[,c('easting','northing')]))
part_mask2 = outer(DF2$year == '2001', DF2$year == '2001') + 
	outer(DF2$year == '2006', DF2$year == '2006')
n2 = dim(D2)[1]
distop2 = D2[1:length(y),(length(y) + 1):n2]
distpp2 = D2[(length(y) + 1):n2,(length(y) + 1):n2]
maskop2 = part_mask2[1:length(y),(length(y) + 1):n2]
maskpp2 = part_mask2[(length(y) + 1):n2,(length(y) + 1):n2]
Xp2 = model.matrix(~ year + dist2road + sideroad:dist2road, data = DFp2)

#strata3
DFp3 = DFpy[DFpy$strat == 3,]
DF3 = rbind(DF[,c('easting','northing','year')],
	DFp3[,c('easting','northing','year')])
D3 = as.matrix(dist(DF3[,c('easting','northing')]))
part_mask3 = outer(DF3$year == '2001', DF3$year == '2001') + 
	outer(DF3$year == '2006', DF3$year == '2006')
n3 = dim(D3)[1]
distop3 = D3[1:length(y),(length(y) + 1):n3]
distpp3 = D3[(length(y) + 1):n3,(length(y) + 1):n3]
maskop3 = part_mask3[1:length(y),(length(y) + 1):n3]
maskpp3 = part_mask3[(length(y) + 1):n3,(length(y) + 1):n3]
Xp3 = model.matrix(~ year + dist2road + sideroad:dist2road, data = DFp3)

maxvar = 4*var(log(y+1))
maxrange = 4*max(distmat)
maxphi = 10000

source(paste0(path,'estpred.R'))
#undebug(estpred)
epout_gamma = estpred(optout_gamma$par, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, 
	maxphi = maxphi, Z1 = Z1, Z2 = Z2, mask = mask, 
	family = 'gamma', Xp = Xp1, 
	dist_op = distop1, dist_pp = distpp1,
	maskop = maskop1, maskpp = maskpp1)
DFp1$w_pred = as.numeric(epout_gamma$w_pred)
DFp1$w_se = as.numeric(epout_gamma$wpred_se_adj)
DFp1_01 = DFp1[DFp1$year == '2001',]
DFp1_06 = DFp1[DFp1$year == '2006',]

epout_gamma = estpred(optout_gamma$par, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, 
	maxphi = maxphi, Z1 = Z1, Z2 = Z2, mask = mask, 
	family = 'gamma', Xp = Xp2, 
	dist_op = distop2, dist_pp = distpp2,
	maskop = maskop2, maskpp = maskpp2)
DFp2$w_pred = as.numeric(epout_gamma$w_pred)
DFp2$w_se = as.numeric(epout_gamma$wpred_se_adj)
DFp2_01 = DFp2[DFp2$year == '2001',]
DFp2_06 = DFp2[DFp2$year == '2006',]

epout_gamma = estpred(optout_gamma$par, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, 
	maxphi = maxphi, Z1 = Z1, Z2 = Z2, mask = mask, 
	family = 'gamma', Xp = Xp3, 
	dist_op = distop3, dist_pp = distpp3,
	maskop = maskop3, maskpp = maskpp3)
DFp3$w_pred = as.numeric(epout_gamma$w_pred)
DFp3$w_se = as.numeric(epout_gamma$wpred_se_adj)
DFp3_01 = DFp3[DFp3$year == '2001',]
DFp3_06 = DFp3[DFp3$year == '2006',]

file_name = "figures/moss_maps"
pdf(paste0(file_name,'.pdf'), width = 15, height = 13)

nclass = 12
cex_1 = .6
cex_2 = 3
cex_3 = 4
mar1 = c(1,0,1,0)
mar2 = c(0,0,0,0)
xleft = .07
xright = xleft + .2
ybot = .2
ytop = .8
brks_cex = 2.2
leg_right = 0.5
textx = -422.5
texty = 2008
xtextleg = 0.23
ytextleg = 0.875
cextextleg = 3.0

# predictions 

layout(matrix(1:6, nrow = 2, byrow = TRUE), widths = c(4,1.2,4))

allbrks = classIntervals(c(DFp1_01$w_pred, DFp1_06$w_pred, 
	DFp2_01$w_pred, DFp2_06$w_pred, DFp3_01$w_pred, DFp3_06$w_pred), 
	n = nclass, style = 'fisher')
	
par(mar = mar1)
plot(DFp1_01[,c('easting','northing')], type = 'n', xaxt = 'n', yaxt = 'n',
	xlab = '', ylab = '')
f10 = classIntervals(DFp3_01$w_pred, n = nclass, style = 'fixed', 
	fixedBreaks = allbrks$brks)
pal = viridis(length(f10$brks) - 1)
f10Colours = findColours(f10, pal)
points(DFp3_01[,c('easting','northing')], pch = 19, col = f10Colours,
	cex = cex_3)
f10 = classIntervals(DFp2_01$w_pred, n = nclass, style = 'fixed', 
	fixedBreaks = allbrks$brks)
pal = viridis(length(f10$brks) - 1)
f10Colours = findColours(f10, pal)
points(DFp2_01[,c('easting','northing')], pch = 19, col = f10Colours,
	cex = cex_2)
f10 = classIntervals(DFp1_01$w_pred, n = nclass, style = 'fixed', 
	fixedBreaks = allbrks$brks)
pal = viridis(length(f10$brks) - 1)
f10Colours = findColours(f10, pal)
points(DFp1_01[,c('easting','northing')], pch = 19, col = f10Colours, 
	cex = cex_1)
text(textx, texty, '2001', cex = 6)

par(mar = mar2)
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = xleft, ybottom = ybot, xright = xright, ytop = ytop,
  breaks = f10$brks, colors = pal, cex = brks_cex, 
  printFormat = "1.2")
text(xtextleg, ytextleg, 'Predictions', cex = cextextleg)

par(mar = mar1)
plot(DFp1_01[,c('easting','northing')], type = 'n', xaxt = 'n', yaxt = 'n',
	xlab = '', ylab = '')
f10 = classIntervals(DFp3_06$w_pred, n = nclass, style = 'fixed', 
	fixedBreaks = allbrks$brks)
pal = viridis(length(f10$brks) - 1)
f10Colours = findColours(f10, pal)
points(DFp3_06[,c('easting','northing')], pch = 19, col = f10Colours,
	cex = cex_3)
f10 = classIntervals(DFp2_06$w_pred, n = nclass, style = 'fixed', 
	fixedBreaks = allbrks$brks)
pal = viridis(length(f10$brks) - 1)
f10Colours = findColours(f10, pal)
points(DFp2_06[,c('easting','northing')], pch = 19, col = f10Colours,
	cex = cex_2)
f10 = classIntervals(DFp1_06$w_pred, n = nclass, style = 'fixed', 
	fixedBreaks = allbrks$brks)
pal = viridis(length(f10$brks) - 1)
f10Colours = findColours(f10, pal)
points(DFp1_06[,c('easting','northing')], pch = 19, col = f10Colours, 
	cex = cex_1)
text(textx, texty, '2006', cex = 6)

# standard errors 

allbrks = classIntervals(c(DFp1_01$w_se, DFp1_06$w_se, 
	DFp2_01$w_se, DFp2_06$w_se, DFp3_01$w_se, DFp3_06$w_se), 
	n = nclass, style = 'fisher')

par(mar = mar1)
plot(DFp1_01[,c('easting','northing')], type = 'n', xaxt = 'n', yaxt = 'n',
	xlab = '', ylab = '')
f10 = classIntervals(DFp3_01$w_se, n = nclass, style = 'fixed', 
	fixedBreaks = allbrks$brks)
pal = magma(length(f10$brks) - 1)
f10Colours = findColours(f10, pal)
points(DFp3_01[,c('easting','northing')], pch = 19, col = f10Colours,
	cex = cex_3)
f10 = classIntervals(DFp2_01$w_se, n = nclass, style = 'fixed', 
	fixedBreaks = allbrks$brks)
pal = magma(length(f10$brks) - 1)
f10Colours = findColours(f10, pal)
points(DFp2_01[,c('easting','northing')], pch = 19, col = f10Colours,
	cex = cex_2)
f10 = classIntervals(DFp1_01$w_se, n = nclass, style = 'fixed', 
	fixedBreaks = allbrks$brks)
pal = magma(length(f10$brks) - 1)
f10Colours = findColours(f10, pal)
points(DFp1_01[,c('easting','northing')], pch = 19, col = f10Colours, 
	cex = cex_1)
points(DF[DF$year == '2001',c('easting','northing')], pch = 4, cex = 6,
	lwd = 2, col = 'green4')
text(textx, texty, '2001', cex = 6)

par(mar = mar2)
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = xleft, ybottom = ybot, xright = xright, ytop = ytop,
  breaks = f10$brks, colors = pal, cex = brks_cex, 
  printFormat = "1.2")
text(xtextleg, ytextleg + .06, 'Standard', cex = cextextleg)
text(xtextleg, ytextleg, 'Errors', cex = cextextleg)

plot(DFp1_01[,c('easting','northing')], type = 'n', xaxt = 'n', yaxt = 'n',
	xlab = '', ylab = '')
f10 = classIntervals(DFp3_06$w_se, n = nclass, style = 'fixed', 
	fixedBreaks = allbrks$brks)
pal = magma(length(f10$brks) - 1)
f10Colours = findColours(f10, pal)
points(DFp3_06[,c('easting','northing')], pch = 19, col = f10Colours,
	cex = cex_3)
f10 = classIntervals(DFp2_06$w_se, n = nclass, style = 'fixed', 
	fixedBreaks = allbrks$brks)
pal = magma(length(f10$brks) - 1)
f10Colours = findColours(f10, pal)
points(DFp2_06[,c('easting','northing')], pch = 19, col = f10Colours,
	cex = cex_2)
f10 = classIntervals(DFp1_06$w_se, n = nclass, style = 'fixed', 
	fixedBreaks = allbrks$brks)
pal = magma(length(f10$brks) - 1)
f10Colours = findColours(f10, pal)
points(DFp1_06[,c('easting','northing')], pch = 19, col = f10Colours, 
	cex = cex_1)
points(DF[DF$year == '2006',c('easting','northing')], pch = 4, cex = 6,
	lwd = 2, col = 'green4')
text(textx, texty, '2006', cex = 6)

layout(1)

dev.off()

system(paste0('pdfcrop ','\'',
	path,file_name,'.pdf','\''))
system(paste0('cp ','\'',
	path,file_name,'-crop.pdf','\' ','\'',
	path,file_name,'.pdf','\''))
system(paste0('rm ','\'',
		path,file_name,'-crop.pdf','\''))
