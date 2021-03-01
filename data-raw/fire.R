library(data.table)
library(magrittr)
library(parallel)
library(plot3D)
options(mc.cores=detectCores())
set.seed(10101)


# get estimated mean charcoal ages of samples in bear and knowles
sid <- char_pmfs %>% colnames
br_kn <- sort(c(grep('BC', sid), grep('LK', sid), grep('UK', sid), grep('DFK', sid), grep('Dam', sid)))
pmfs <- char_pmfs %>% as.matrix
pmfs <- pmfs[, br_kn]
char <- charcoal[br_kn, ]
br_ca <- max(creeks$contr_area[creeks$creek_name == 'bear'])
kn_ca <- max(creeks$contr_area[creeks$creek_name == 'knowles'])
# df <- pmfs[, char$facies == 'DF']

# order by weighted mean
index <- char_pmfs %>% rownames %>% as.numeric %>% sort
dfmn <- apply(df, 2, function (x) weighted.mean(rev(index), x))
dfmns <- apply(dfs, 1, function (x) weighted.mean(rev(index), x))
charmn <- apply(pmfs, 2, function (x) weighted.mean(rev(index), x))



# stream input/output rate - global
tis <- c(.01, .20)
tos <- c(.01, .20)

begin <- Sys.time()
diff <- as.difftime('08:00:00', '%H:%M:%S', units = 'hour')
(end <- begin + diff)
i <- 0
while (Sys.time() < end) {
  res <- accumulater(nrow(dfs), tis, tos, ia_ar, dfmns, batch = 100, it = 500)
  rec <- rbind(rec, res)
  save(rec, file = "fire_io_20201027.rds")
}

load('fire_io_20201027.rds')
# rec <- res

(good <- rec[rec$ks == min(rec$ks), ])

# png('ia_top.png', height = 20, width = 20, units = 'cm', res = 300)
scatter3D(rec[,1], rec[,2], rec[,4], ticktype = 'detailed', pch = 20,
          xlab = "input rate", ylab = "output rate", zlab = "K-S stat",
          phi = 30, theta = 120)
# dev.off()

kpmin <- rec[rec$kp == min(rec$kp), ]
ksmin <- rec[rec$ks == min(rec$ks), ]
png('fit_ia.png', height = 14, width = 17, units = 'cm', res = 300)
plot(emp_cdf(dfmns), xlim = c(0, 11000),
     xlab = 'charcoal age', ylab = 'CDF')
for (i in 1:nrow(ksmin)) {
  lines(emp_cdf(accumulate(nrow(dfs), ksmin[i,1], ksmin[i,2], ia_ar, 1000)),
        lwd = 2, col = get_palette('crimson', .33))
}
for (i in 1:nrow(kpmin)) {
  lines(emp_cdf(accumulate(nrow(dfs), kpmin[i,1], kpmin[i,2], ia_ar, 1000)),
        lwd = 2, col = get_palette('ocean', .33))
}
legend('bottomright', legend = c('observed', 'K-S fit', 'Kuiper fit'),
       fill = get_palette(c('charcoal', 'crimson', 'ocean'), .7))
dev.off()

# charcoal input rate - global
# tis <- c(.01, .25)
# tos <- c(.01, .25)
# tis <- c(.14, .19)
# tos <- c(.12, .17)

begin <- Sys.time()
diff <- as.difftime('01:00:00', '%H:%M:%S', units = 'hour')
(end <- begin + diff)
i <- 0
while (Sys.time() < end) {
  res <- accumulater(length(charmn), tis, tos, 0, charmn, 100, 500)
  rec <- rbind(rec, res)
  save(rec, file = "fire_20201027.rds")
}

# load('fire_20201027.rds')
# rec <- res

(good <- rec[rec$ks == min(rec$ks), ])

# png('fire_ti.png', height = 17, width = 15, units = 'cm', res = 300)
scatter3D(rec[,1], rec[,2], rec[,3], ticktype = 'detailed', pch = 20,
          xlab = "input rate", ylab = "output rate", zlab = "test value",
          phi = 30, theta = 120)
# dev.off()

glti_mn <- mean(good$ti)
glto_mn <- mean(good$to)


# bootstrap confidence intervals - global
# tis <- c(.01, .25)
# tos <- c(.01, .25)
begin <- Sys.time()
diff <- as.difftime('09:00:00', '%H:%M:%S', units = 'hour')
(end <- begin + diff)
i <- 0
while (Sys.time() < end) {
  res <- boot_accum(length(charmn), tis, tos, 0, charmn, 1000, 100)
  rec <- rbind(rec, res)
  save(rec, file = "fireboot_20201029.rds")
}

# load('fireboot_20201027.rds')
# rec <- res

(good <- rec[rec$ks == min(rec$ks), ])
nrow(rec)

scatter3D(rec[,1], rec[,2], rec[,3], ticktype = 'detailed', pch = 20,
          xlab = "input rate", ylab = "output rate", zlab = "test value",
          phi = 30, theta = 145)

load('fireboot_20201029.rds')
gl_rates <- data.frame(rate = rec$ti / (br_ca + kn_ca), type = 'input', creek = 'global')
gl_rates <- rbind(gl_rates, data.frame(
  rate = rec$to / (br_ca + kn_ca), type = 'output', creek = 'global'
))

boxplot(rate ~ type, data = gl_rates)


gl_boot <- rec[order(rec$ti), ]
glti_lwr <- gl_boot[floor(nrow(rec) * .025), ]
glti_upr <- gl_boot[ceiling(nrow(rec) * .975), ]

gl_boot <- rec[order(rec$to), ]
glto_lwr <- gl_boot[floor(nrow(rec) * .025), ]
glto_upr <- gl_boot[ceiling(nrow(rec) * .975), ]



tis <- c(.001, .25)
tos <- c(.05, .13)

# get estimated mean charcoal ages of samples in bear and knowles
sid <- char_pmfs %>% colnames
br <- sort(c(grep('BC', sid)))
kn <- sort(c(grep('LK', sid), grep('UK', sid), grep('DFK', sid), grep('Dam', sid)))
pmfs <- char_pmfs %>% as.matrix
brpmfs <- pmfs[, br]
knpmfs <- pmfs[, kn]
brchar <- charcoal[br, ]
knchar <- charcoal[kn, ]
brdf <- brpmfs[, brchar$facies == 'DF']
kndf <- knpmfs[, knchar$facies == 'DF']

# order by weighted mean
index <- char_pmfs %>% rownames %>% as.numeric %>% sort
brdfmn <- apply(brdf, 2, function (x) weighted.mean(rev(index), x))
kndfmn <- apply(kndf, 2, function (x) weighted.mean(rev(index), x))
brdfmns <- apply(brpmfs, 2, function (x) weighted.mean(rev(index), x))
kndfmns <- apply(knpmfs, 2, function (x) weighted.mean(rev(index), x))

tis <- c(.001, .25)
tos <- c(.001, .25)

# estimate input rate of charcoal from fire for bear

# stream input/output rate corrected for inherited age
# bear has no inherited age samples, using inherited ages from knowles
tis <- c(.001, .2)
tos <- c(.001, .2)

begin <- Sys.time()
diff <- as.difftime('02:00:00', '%H:%M:%S', units = 'hour')
(end <- begin + diff)
i <- 0
while (Sys.time() < end) {
  res <- accumulater(ncol(brpmfs), tis, tos, ia_ar, brdfmns, 100, 500)
  rec <- rbind(rec, res)
  save(rec, file = "fire_brio_20201030.rds")
}

# load('fire_brio_20201030.rds')
# rec <- res

(good <- rec[rec$ks == min(rec$ks), ])

scatter3D(rec[,1], rec[,2], rec[,3], ticktype = 'detailed', pch = 20,
          xlab = "input rate", ylab = "output rate", zlab = "test value",
          phi = 30, theta = 320)


brti_mn <- mean(good$ti)
brto_mn <- mean(good$to)



# bootstrap confidence intervals for bear
begin <- Sys.time()
diff <- as.difftime('04:00:00', '%H:%M:%S', units = 'hour')
(end <- begin + diff)
i <- 0
while (Sys.time() < end) {
  res <- boot_accum(ncol(brpmfs), tis, tos, 0, brdfmns, 1000, 100)
  rec <- rbind(rec, res)
  save(rec, file = "fireboot_br_20201030.rds")
}

# load('fireboot_br_20201030.rds')
# rec <- res

(good <- rec[rec$ks == min(rec$ks), ])
nrow(rec)

scatter3D(rec[,1], rec[,2], rec[,3], ticktype = 'detailed', pch = 20,
          xlab = "input rate", ylab = "output rate", zlab = "test value",
          phi = 30, theta = 145)

br_boot <- rec[order(rec$ti), ]
brti_lwr <- br_boot[floor(nrow(rec) * .025), ]
brti_upr <- br_boot[ceiling(nrow(rec) * .975), ]

br_boot <- rec[order(rec$to), ]
brto_lwr <- br_boot[floor(nrow(rec) * .025), ]
brto_upr <- br_boot[ceiling(nrow(rec) * .975), ]

br_rates <- data.frame(rate = rec$ti / br_ca, type = 'input', creek = 'bear')
br_rates <- rbind(br_rates, data.frame(
  rate = rec$to / br_ca, type = 'output', creek = 'bear'
))

boxplot(rate ~ type, data = br_rates)

rates <- rbind(gl_rates, br_rates)
boxplot(rate ~ type + creek, data = rates)



# estimate input rate of charcoal from fire for knowles
tis <- c(.001, .25)
tos <- c(.001, .25)

begin <- Sys.time()
diff <- as.difftime('02:00:00', '%H:%M:%S', units = 'hour')
(end <- begin + diff)
i <- 0
while (Sys.time() < end) {
  res <- accumulater(length(kndfmns), tis, tos, 0, kndfmns, 100, 500)
  rec <- rbind(rec, res)
  save(rec, file = "fire_knio_20201101.rds")
}

# load('fire_knio_20201024.rds')
# rec <- res
nrow(rec)

(good <- rec[rec$ks == min(rec$ks), ])

scatter3D(rec[,1], rec[,2], rec[,3], ticktype = 'detailed', pch = 20,
          xlab = "input rate", ylab = "output rate", zlab = "test value",
          phi = 30, theta = 320)

knti_mn <- mean(good$ti)
knto_mn <- mean(good$to)


setwd('/home/crumplecup/work/')
# bootstrap confidence intervals for knowles
begin <- Sys.time()
diff <- as.difftime('16:00:00', '%H:%M:%S', units = 'hour')
(end <- begin + diff)
i <- 0
while (Sys.time() < end) {
  res <- boot_accum(ncol(knpmfs), tis, tos, 0, kndfmns, 1000, 100)
  rec <- rbind(rec, res)
  save(rec, file = "fireboot_kn_20201031.rds")
}

# load('fireboot_kn_20201031.rds')
# rec <- res
nrow(rec)

(good <- rec[rec$ks == min(rec$ks), ])

scatter3D(rec[,1], rec[,2], rec[,3], ticktype = 'detailed', pch = 20,
          xlab = "input rate", ylab = "output rate", zlab = "test value",
          phi = 30, theta = 145)

kn_boot <- rec[order(rec$ti), ]
knti_lwr <- kn_boot[floor(nrow(rec) * .025), ]
knti_upr <- kn_boot[ceiling(nrow(rec) * .975), ]

kn_boot <- rec[order(rec$to), ]
knto_lwr <- kn_boot[floor(nrow(rec) * .025), ]
knto_upr <- kn_boot[ceiling(nrow(rec) * .975), ]

kn_rates <- data.frame(rate = rec$ti / kn_ca, type = 'input', creek = 'knowles')
kn_rates <- rbind(kn_rates, data.frame(
  rate = rec$to / kn_ca, type = 'output', creek = 'knowles'
))

boxplot(rate ~ type, data = kn_rates)

pal <- get_palette(c('ocean', 'hardwood', 'forest'))
rates <- rbind(gl_rates, br_rates, kn_rates)
png('creeks_io.png', height = 15, width = 17, units = 'cm', res = 300)
boxplot(rate ~ type + creek, data = rates, axes = F,
        xlab = 'creek', ylab = 'rate per km2 CA',
        col = pal[c(1,1,2,2,3,3)])
axis(side = 2)
axis(side = 1,
     at = c(1:6), labels = c(
       "in", 'out',
       "in", 'out',
       'in', 'out'))
legend("bottomright", legend = c('global', 'bear', 'knowles'),
       fill = pal, horiz = F)
dev.off()

knti_lwr / kn_ca
brti_upr / br_ca
brti_lwr / br_ca


# convert rate per study area to rate to km2

# number of nodes per study area
br_N <- nrow(creeks[creeks$creek_name == 'bear', ])
kn_N <- nrow(creeks[creeks$creek_name == 'knowles', ])

# number of sample sites per study area
ct <- charcoal[ , .N, by = c('facies', 'family')]
br_n <- ct[grep('BC', ct$family), ] %>% nrow
kn_n <- ct[
  c(grep('UK', ct$family),
  grep('LK', ct$family),
  grep('DFK', ct$family),
  grep('Dam', ct$family)),
  ] %>% nrow


# 3D plot of delivery weighting function

crks <- creeks[creeks$creek_name %in% c('bear', 'knowles'), ]
crks <- rater1(crks, .11, .09)
crks$elev <- slope_to_elev(crks$slope, crks$ToMouth_km)
crks$elev[crks$creek_name == 'knowles'] <- crks$elev[crks$creek_name == 'knowles'] - min(crks$elev[crks$creek_name == 'knowles'])
crk_xy <- coordinates(crks)

crks$wt_cols <- 0
crks$wt_cols[crks$dp <= 0.007] <- get_palette('ocean', .1)
crks$wt_cols[crks$dp > 0.007 & crks$dp <= .01] <- get_palette('sky', .1)
crks$wt_cols[crks$dp > .01 & crks$dp <= .015] <- get_palette('gold', .1)
crks$wt_cols[crks$dp > .015 & crks$dp <= .02] <- get_palette('rose', .1)
crks$wt_cols[crks$dp > .02] <- get_palette('crimson', .1)

png('crk_wts.png', height = 17, width = 17, units = 'cm', res = 300)
scatter3D(crk_xy[,1], crk_xy[,2], crks$elev, col = crks$wt_cols, pch= 20,
          phi = 20, theta = 150)
legend('bottomright', legend = c('wt < 1', '1 < wt < 2', 'wt > 2'),
       fill = get_palette(c('ocean', 'gold', 'crimson'), .7))
dev.off()

library(rgdal)
writeOGR(crks, '/media/crumplecup/catacomb/gis/creeks/kn_br_wts.shp', 'creeks', 'ESRI Shapefile')


load('/media/crumplecup/catacomb/work/vol_20201008.rds')

(good <- rec[rec$ks_lvl + rec$ks_vol == min(rec$ks_lvl + rec$ks_vol), ])
best_fit <- as.numeric(good)
kn <- creeks[creeks$creek_name == 'knowles', ]
best_so <- rater1(kn, best_fit[1], best_fit[2])

pred <- backfill1(best_so, 10000, best_fit[3], best_fit[4], best_fit[5], best_fit[6])
# pred <- fit_volumes1(best_so, 10000, best_fit[3], best_fit[4], 10)
pred[[2]]
pred <- pred[[1]]

# setwd('/home/crumplecup/work')
# png('vol_fit_10-08.png', height = 15, width = 17, units = 'cm', res = 300)
plot(kn$ToMouth_km, kn$xsec_area, type = 'l', lwd = 2.5, col = get_palette('forest', .7),
     xlab = 'distance to outlet (km)', ylab = 'unit cross-sectional area (km2)')
lines(kn$ToMouth_km, rev(pred$vol), lwd = 2, col = get_palette('charcoal', .7))
legend('topright', legend = c('observed', 'fit'),
       fill = get_palette(c('forest', 'charcoal'), .7))
# dev.off()

# png('xsec_fit_10-08.png', height = 15, width = 17, units = 'cm', res = 300)
plot(best_so$xsec_area %>% emp_cdf, pch = 20, col = get_palette('ocean', .2),
     xlab = 'unit cross-sectional area m2', ylab = 'CDF')
lines(pred$vol %>% emp_cdf, lwd = 3, col = get_palette('charcoal', .6))
legend('bottomright', legend = c('observed', 'fit'),
       fill = get_palette(c('ocean', 'charcoal'), .7))
# dev.off()

# png('lvl_fit_10-08.png', height = 15, width = 17, units = 'cm', res = 300)
plot(pred$lvl %>% emp_cdf, type = 'l', lwd = 3, col = get_palette('charcoal', .7),
     xlab = 'average valley depth m2', ylab = 'CDF', xlim = c(0, max(best_so$lvl)))
points(best_so$lvl %>% emp_cdf, pch = 20, col = get_palette('ocean', .2))
legend('bottomright', legend = c('observed', 'fit'),
       fill = get_palette(c('ocean', 'charcoal'), .7))
# dev.off()

# CDF of exp distr is 1 - e^lambda*t

old <- 5000
gli_cdf <- 1 - exp(-glti_mn * 0:old)
kni_cdf <- 1 - exp(-knti_mn * 0:old)
bri_cdf <- 1 - exp(-brti_mn * 0:old)
glo_cdf <- 1 - exp(-glto_mn * 0:old)
kno_cdf <- 1 - exp(-knto_mn * 0:old)
bro_cdf <- 1 - exp(-brto_mn * 0:old)

kn_it <- 0
for (i in seq_along(kndfmns)) {
  if (i > 1) {
    kn_it[i-1] <- kndfmns[i] - kndfmns[i-1]
  }
}

br_it <- 0
for (i in seq_along(brdfmns)) {
  if (i > 1) {
    br_it[i-1] <- brdfmns[i] - brdfmns[i-1]
  }
}

glmns <- sort(c(kndfmns, brdfmns))
gl_it <- 0
for (i in seq_along(glmns)) {
  if (i > 1) {
    gl_it[i-1] <- glmns[i] - glmns[i-1]
  }
}


# png('char_cdfs.png', height = 17, width = 23, units = 'cm', res = 300)
plot(gli_cdf, log = 'x', ylim = c(0,1), xlim = c(1, 1000),
     type = 'l', lwd = 2, col = pal[1],
     xlab = 't', ylab = 'CDF')
lines(glo_cdf, lwd = 2, lty = 2, col = pal[1])
lines(bri_cdf, lwd = 2, col = pal[2])
lines(bro_cdf, lwd = 2, lty = 2, col = pal[2])
lines(kni_cdf, lwd = 2, col = pal[3])
lines(kno_cdf, lwd = 2, lty = 2, col = pal[3])
points(emp_cdf(gl_it), pch = 20, col = pal[1])
points(emp_cdf(br_it), pch = 20, col = pal[2])
points(emp_cdf(kn_it), pch = 20, col = pal[3])
legend('bottomright', legend = c('global', 'bear', 'knowles', 'in', 'out', 'observed'),
       fill = c(pal, NA, NA, NA), border = c(NA),
       lty = c(NA, NA, NA, 1, 2, NA),
       pch = c(NA, NA, NA, NA, NA, 1))
# dev.off()




# back of the envelope lambda investigation

t <- 1:100
lin <- 0.11
lout <- 0.08
k <- 1:10

poisson_pmf <- function(l, t, k) {
  ((l * t)^k * exp(-l * t)) / factorial(k)
}

dual_poisson <- function(lin, lout, t, k) {
  ((lin * (1 - lout) * t)^k - exp(-(lin * 1 - lout) * t)) / factorial(k)
}

dual_exp <- function(lin, lout, t) {
  lin * t * (1 - lout) * t * exp(-(lin * (1 - lout) * t * t))
}

dual_cdf <- function(lin, lout, t) {
  1 - exp(-(lin * (1 - lout) * t * t))
}

exp_cdf <- function(l, t) {
  1 - exp(-(l * t))
}

bin <- poisson_pmf(lin, t, 1)
bout <- poisson_pmf(lout, t, 1)

pal <- get_palette(c('ocean', 'forest', 'gold', 'crimson'))

plot(bout, col = pal[1])
lines(bin, lwd = 3, col = pal[2])
lines(bit, lwd = 3, col = pal[3])
plot(poisson_pmf(lin, t, 1) - poisson_pmf(lout, t, 1))
plot(dual_poisson(lin, lout, t, 1))
plot(dual_exp(lin, lout, t))
plot(dual_cdf(lin, lout, t))
plot(exp_cdf(lin, t))
plot(exp_cdf(.5, 1:10))
plot(exp_cdf(.25, 1:20))

library(extraDistr)

plot(density(rskellam(10000, .11, .08)))
plot(dskellam(-3:3, .11, .08))
dskellam(-3:3, .11, .08)
dskellam(-3:3, .011, .008)
bit <- dskellam(0:100, 110, 80)
plot(cumsum(bit))

bot <- skellam::dskellam(-20:20, 10, 5)
plot(-20:20, bot)
