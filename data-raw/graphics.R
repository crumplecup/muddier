
# subset tyee formation from oregon geologic map
library(rgdal)

fgdb <- '/media/crumplecup/Seagate Backup Plus Drive/gis/or_geology/OGDC_v6.gdb'

subset(ogrDrivers(), grepl("GDB", name))
fc_list <- ogrListLayers(fgdb)
print(fc_list)

fc <- readOGR(dsn=fgdb,layer = fc_list[1])
names(fc)
tyee <- fc[grep('Tyee', fc$FORMATION), ]
tyee <- rgeos::gUnaryUnion(tyee)

tyee <- sp::SpatialPolygonsDataFrame(tyee, data.frame(name = 'tyee'), F)

rgdal::writeOGR(tyee,
                '/media/crumplecup/Seagate Backup Plus Drive/gis/or_geology/',
                'tyee',
                driver = 'ESRI Shapefile')

#subset rivers of the study area from NHD HR set
gis_dir <- '/media/crumplecup/Seagate Backup Plus Drive/gis'
nhd_dir <- file.path(gis_dir, 'nhd_hr/NHDPLUS_H_1710_HU4_GDB')
nhd_gdb <- file.path(nhd_dir, 'NHDPLUS_H_1710_HU4_GDB.gdb')
nhd_list <- ogrListLayers(nhd_gdb)
print(nhd_list)

nhd <- readOGR(dsn = nhd_gdb, layer = nhd_list[56])
names(nhd)
sius <- nhd[grep('Siuslaw', nhd$GNIS_Name), ]
ump <- nhd[grep('Umpqua', nhd$GNIS_Name), ]
kno <- nhd[grep('Knowles', nhd$GNIS_Name), ]
hof <- nhd[grep('Hoffman', nhd$GNIS_Name), ]
gol <- nhd[grep('Golden', nhd$GNIS_Name), ]
hoff <- nhd[grep(17100206000368, nhd$ReachCode), ] # trib of Hoffman, study area

ced <- nhd[nhd$GNIS_ID == '01139500', ] # Cedar Creek
swe <- nhd[nhd$GNIS_ID == '01150779', ] # Sweet Creek
kno <- nhd[nhd$GNIS_ID == '01144693', ] # Knowles Creek
smi <- nhd[nhd$GNIS_ID == '01149763', ] # Smith Creek
was <- nhd[nhd$GNIS_ID == '01639342', ] # Wasson Creek
siu <- nhd[nhd$GNIS_ID == '01149557', ] # Siuslaw River

bea <- nhd[nhd$ReachCode %in% c(17100206000419), ] # Bear Creek

ceda <- nhd[nhd$ReachCode %in% c(
  17100206009827,
  17100206009832,
  17100206009834,
  17100206009837,
  17100206009851,
  17100206035016
), ] # Cedar Creek trib, study area

gol <- nhd[nhd$ReachCode %in% c(
  17100303000683,
  17100303000698,
  17100303000709,
  17100303000730,
  17100303004325,
  17100303034656
), ] # Golden Ridge Creek

gold <- nhd[nhd$ReachCode %in% c(
  17100303000699
), ] # Golden Ridge Creek trib

golde <- nhd[nhd$ReachCode %in% c(
  17100303003877,
  17100303034621
), ] # Golden Ridge Creek trib

work_dir <- '/home/crumplecup/work'
writeOGR(siu, work_dir, 'siuslaw river', driver = 'ESRI Shapefile', overwrite = T)
writeOGR(hof, work_dir, 'hoffman creek', driver = 'ESRI Shapefile', overwrite = T)
writeOGR(hoff, work_dir, 'hoffman trib', driver = 'ESRI Shapefile', overwrite = T)
writeOGR(kno, work_dir, 'knowles creek', driver = 'ESRI Shapefile', overwrite = T)
writeOGR(bea, work_dir, 'bear creek', driver = 'ESRI Shapefile', overwrite = T)
writeOGR(swe, work_dir, 'sweet creek', driver = 'ESRI Shapefile', overwrite = T)
writeOGR(ced, work_dir, 'cedar creek', driver = 'ESRI Shapefile', overwrite = T)
writeOGR(ceda, work_dir, 'cedar trib', driver = 'ESRI Shapefile', overwrite = T)
writeOGR(ump, work_dir, 'umpqua river', driver = 'ESRI Shapefile', overwrite = T)
writeOGR(smi, work_dir, 'smith river', driver = 'ESRI Shapefile', overwrite = T)
writeOGR(was, work_dir, 'wasson creek', driver = 'ESRI Shapefile', overwrite = T)
writeOGR(gol, work_dir, 'golden ridge creek', driver = 'ESRI Shapefile', overwrite = T)
writeOGR(gold, work_dir, 'grc trib1', driver = 'ESRI Shapefile', overwrite = T)
writeOGR(golde, work_dir, 'grc trib2', driver = 'ESRI Shapefile', overwrite = T)


know <- creeks[creeks$creek_name %in% 'knowles', ]
names(nhd)
knowl <- nhd[nhd$Permanent_Identifier %in% c(
  161949916,
  161949915,
  161949899,
  161949849,
  161950029,
  161950028,
  161950070,
  161950069,
  161949912,
  161949950,
  161949949,
  161949995,
  161949994,
  161949993,
  161949906,
  161949905,
  161949814,
  161950013,
  161950012,
  161949904,
  161949792,
  161950041,
  161950040,
  161949776,
  161950054
), ]

knowl <- Line(coordinates(know))
knowl <- Lines(list(knowl), ID = 'knowles')
knowl <- SpatialLines(list(knowl), crs_ref)
knowl <- SpatialLinesDataFrame(knowl, data.frame(ID = 'knowles'), match.ID = F)
knowl <- sp::spTransform(knowl, raster::crs(siu))
writeOGR(knowl, work_dir, 'knowles_study_area', driver = 'ESRI Shapefile', overwrite = T)

plot(knowl)
plot(sp::spTransform(tyee, crs_ref),
     col = get_palette('coral', .1), border = get_palette('coral'))
lines(sp::spTransform(siu, crs_ref), col = get_palette('ocean', .1))
lines(sp::spTransform(hof, crs_ref), col = get_palette('ocean', .1))
lines(sp::spTransform(kno, crs_ref), col = get_palette('ocean', .1))
lines(sp::spTransform(ced, crs_ref), col = get_palette('ocean', .1))
lines(sp::spTransform(swe, crs_ref), col = get_palette('ocean', .1))
lines(sp::spTransform(ump, crs_ref), col = get_palette('ocean', .1))
lines(sp::spTransform(smi, crs_ref), col = get_palette('ocean', .1))
lines(sp::spTransform(was, crs_ref), col = get_palette('ocean', .1))
lines(sp::spTransform(gol, crs_ref), col = get_palette('ocean', .1))
lines(sp::spTransform(hoff, crs_ref), col = get_palette('forest', .1))
lines(sp::spTransform(bea, crs_ref), col = get_palette('forest', .1))
lines(sp::spTransform(ceda, crs_ref), col = get_palette('forest', .1))
lines(sp::spTransform(gold, crs_ref), col = get_palette('forest', .1))
lines(sp::spTransform(golde, crs_ref), col = get_palette('forest', .1))


plot(gol)
plot(was, add = T)
plot(gold, add = T)
plot(golde, add = T)


plot(was)
plot(smi, add = T)
plot(sp::spTransform(was, raster::crs(creeks)))
plot(creeks[creeks$creek_name == 'grc', ], add = T)

plot(swe)
plot(ced, add = T)
plot(ceda, add = T)
plot(kno, add = T)
swe <- nhd[grep(17100206000019, nhd$ReachCode), ] # Sweet Creek

cd1 <- nhd[grep(17100206009827, nhd$ReachCode), ] # cedar trib frag

plot(c(landmarks$x, landmarks$y), add = T)
plot(sp::spTransform(kno, raster::crs(creeks)))
plot(creeks[creeks$creek_name %in% c('bear', 'knowles'), ], add = T)


mydir <- '/media/crumplecup/catacomb/gis/benton_gis'
setwd(mydir)
home <- readOGR('zoning_draft.shp')
head(home, 100)

city_limits <- home[grep('City', home$LAYER), ]
writeOGR(city_limits, file.path(mydir, 'city_limits.shp'), 'city_limits',
                                driver = 'ESRI Shapefile')
ugb <- readOGR('philomath_ugb_draft.shp')
plot(ugb)
lines(city_limits, col = get_palette('crimson'))

plot(hof)
plot(hoff, add = T)
plot(sius, add = T)
plot(sius)
plot(swe, add = T)
plot(swe)
plot(ced, add = T)
plot(ced)

y <- c(0, 1, 2, 4, 8, 16, 32, 64, 96, 112, 120, 124, 126, 127,
       126, 124, 120, 112, 96, 64, 48, 40, 36, 34, 33)
y <- c(y, rev(y))
y <- y / sum(y)
x <- 1:length(y) + 555
setwd(work_dir)
png('eg_tab_1.png')
plot(x, y, type = 'l', lwd = 3, col = get_palette('ocean', .6),
     xlab = 'sample age', ylab = 'probability')
dev.off()

ys <- c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096)
ys <- c(ys, cumsum(rev(ys))[-1])
ys <- c(ys, rev(ys))
ys <- ys / sum(ys)
xs <- 1:length(ys) +175
plot(xs, ys)

eg_tab <- matrix(c(x, y, xs, ys), nrow = 4, byrow = T)
write.csv(eg_tab, 'eg_tab.csv')

png('eg_tab_2.png')
plot(xs, ys, type = 'l', lwd = 3, col = get_palette('ocean', .6),
     xlab = 'sample age', ylab = 'probability')
dev.off()

png('eg_tab_3.png')
plot(x, y, type = 'l', lwd = 3, col = get_palette('ocean', .6),
     xlab = 'sample age', ylab = 'probability')
points(x, y, col = get_palette('charcoal', 1))
dev.off()

png('eg_tab_4.png')
plot(xs, ys, type = 'l', lwd = 3, col = get_palette('ocean', .6),
     xlab = 'sample age', ylab = 'probability')
points(xs, ys, col = get_palette('charcoal', 1))
dev.off()

xar <- array(0, c(length(x),length(xs)))
yar <- xar
xv <- vector(length = length(x)^2, mode = 'numeric')
yv <- vector(length = length(y)^2, mode = 'numeric')

k <- 1
for(i in 1:length(x)) {
  for(j in 1:length(xs)) {
    xv[k] <- x[i] - xs[j]
    yv[k] <- y[i] * ys[j]
    k <- k + 1
  }
}

sum(yv)

?convo
index <- 1:1000
yv <- vector(length = 1000, mode = 'numeric')
ysv <- vector(length = 1000, mode = 'numeric')
for (i in seq_along(x)) {
  yv[x[i]] <- y[i]
  ysv[xs[i]] <- ys[i]
}

bit <- convo(yv, ysv, index)

png('eg_res.png')
plot(bit, xlim = c(345, 415), type = 'l', lwd = 3, col = get_palette('ocean', .6),
     xlab = 'age', ylab = 'probability')
points(bit)
dev.off()


xar <- array(0, c(length(x), length(xs)))
yar <- xar

for (i in seq_along(x)) {
  for (j in seq_along(xs)) {
    xar[i,j] <- x[i] - xs[j]
    yar[i,j] <- y[i] * ys[j]
  }
}

xar <- as.data.frame(xar)
yar <- as.data.frame(yar)

rownames(xar) <- x
colnames(xar) <- xs
rownames(yar) <- x
colnames(yar) <- xs

write.csv(xar, file = 'eg_xar.csv')
write.csv(yar, file = 'eg_yar.csv')

data('df_ph_cdf', package = 'muddier')
data('ff_ph_cdf', package = 'muddier')
data('fg_ph_cdf', package = 'muddier')

pal <- get_palette(c('sky', 'ocean', 'gold', 'coral', 'charcoal', 'forest'), .7)
png('fit_inherited_age.png', height = 6, width = 6, units = 'in', res = 300)
plot(index, df_cdf, type = 'l', lwd = 3, col = pal[5], xlim = c(0, 10000),
     xlab = 'Inherited Age', ylab = 'proportion at or above age(x)')
lines(index, df_ph_cdf, lwd = 3, col = pal[6], lty = 2)
lines(index, fg_cdf, lwd = 3, col = pal[3])
lines(index, fg_ph_cdf, lwd = 3, col = pal[4], lty = 2)
lines(index, ff_cdf, lwd = 3, col = pal[1])
lines(index, ff_ph_cdf, lwd = 3, col = pal[2], lty = 2)
legend('bottomright', legend = c(
  'obs FF', 'fit FF', 'obs FG', 'fit FG', 'obs DF', 'fit DF'),
  fill = pal)
dev.off()

