library(spData)
library(sf)

data(elect80)

DF = elect80[substr(elect80$FIPS, start = 1, stop = 2) == '48',]
dim(DF)[1]

DF = st_as_sf(DF, 4326)
DF = st_transform(DF, crs = 5070)
xy = st_coordinates(DF)/1000
DF$bin_turnout = (DF$pc_turnout > 0.5)*1
Tex_turnout80 = DF
path = '/mnt/ExtraDrive1/Work/desktop_data/2022_papers/hglmm/'
save(Tex_turnout80, file = paste0(path,'Tex_turnout80.rda'))

layout(matrix(1:2, nrow = 1), width = c(5,1))
f10 = classIntervals(logit(DF$pc_turnout), n = 10, style = 'fisher')
pal = viridis(length(f10$brks) - 1)
f10Colours = findColours(f10, pal)
plot(st_geometry(DF), col = f10Colours, pch = 19, cex = 2)
text(-961569, 1503770, 'D', cex = lab_cex, adj = c(0,1))
par(mar = c(0,0,0,0))
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .1, xright = .2, ytop = .9,
  breaks = f10$brks, colors = pal, cex = brks_cex, 
  printFormat = "1.2")

