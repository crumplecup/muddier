library(magrittr)
setwd('/media/erik/catacomb/research')

rec <- data.table::fread('/home/erik/output/df20k_1000.csv')
rec1 <- data.table::fread('/home/erik/output/df20k_1001.csv')
recf <- data.table::fread('/home/erik/output/ff20k_1000.csv')
recf1 <- data.table::fread('/home/erik/output/ff20k_1001.csv')
recf2 <- data.table::fread('/home/erik/output/ff20k_1002.csv')
recg <- data.table::fread('/home/erik/output/fg20k_1000.csv')
recg1 <- data.table::fread('/home/erik/output/fg20k_1001.csv')
rec <- rbind(rec, rec1)
recf <- rbind(recf, rbind(recf1, recf2))
recg <- rbind(recg, recg1)
rm(rec1, recf1, recf2, recg1)

rec10 <- data.table::fread('/home/erik/output/steady_10ky_20kr_7777.csv')

plot(rec$input, rec$n, pch = 20, col = get_palette('charcoal'))
points(rec10$input, rec10$n, pch = 20, col = get_palette('crimson'))


df_stereo <- data.table::fread('/home/erik/output/df_stereo.csv') %>% unlist

# knowles creek pebble counts
peb <- data.table::fread('/media/erik/catacomb/research/knowles_pebble_count.csv')
pebs <- data.table::fread('/media/erik/catacomb/research/knowles_subsurface_pebbles.csv')

# debris-flow transit times at selected fits
df_transit_ad <- data.table::fread('/home/erik/output/df_transits_ad.csv') %>% unlist
df_transit_ch <- data.table::fread('/home/erik/output/df_transits_ch.csv') %>% unlist
df_transit_kp <- data.table::fread('/home/erik/output/df_transits_kp.csv') %>% unlist
df_transit_ks <- data.table::fread('/home/erik/output/df_transits_ks.csv') %>% unlist


cdf_hook <- function(obs) {
  # given observations in random order
  # returns a vector of the cdf value of each ob
  vals <- rep(0, length(obs))
  cdf <- emp_cdf(obs)
  x <- cdf[ , 1]
  y <- cdf[ , 2]
  for (i in 1:length(x)) {
    if (i == 1) {
      vals[obs <= x[i]] <- y[i]
    }
    if (i > 1) {
      vals[obs > x[i-1] & obs <= x[i]] <- y[i]
    }
  }
  vals
}

rarer <- function(d, d50, mod, win = .95, r = 1) {
  # given a b-axis width, the d50, and a model fit
  # returns the relative rarity of the diameter

  # intercept at surface d50
  off <- log(d50) * mod$coefficients[2] - log(0.5)
  # pred width when less common
  upr <- exp((off + log(0.5 + 0.5 * win)) / mod$coefficients[2])
  # pred width when less rare
  lwr <- exp((off + log(0.5 - 0.5 * win)) / mod$coefficients[2])

  # test if d is in range, else recurse up or down
  if (d < lwr) {
    r <- rarer(d, lwr, mod, r = r * 0.5 / (0.5 - 0.5 * win))
  }

  if (d > upr) {
    r <- rarer(d, upr, mod, r = r * 0.5 / (0.5 + 0.5 * win))
  }

  # cdf at d
  c <- exp(log(d) * mod$coefficients[2] - off)
  # relative rarity to d50
  r <- r * 0.5 / c
  r
}




# kolmogorov-smirnov
# debris-flow input/output rate
drt <- rec$input[rec$ks == min(rec$ks)[1]]
drt_ad <- rec$input[rec$ad == min(rec$ad)[1]]
drt_ch <- rec$input[rec$ch == min(rec$ch)[1]]
drt_kp <- rec$input[rec$kp == min(rec$kp)[1]]
drt_ks <- rec$input[rec$ks == min(rec$ks)[1]]
# fines input/output rate
frt <- recf$input[recf$ks == min(recf$ks)[1]][1]
frt_ad <- recf$input[recf$ad == min(recf$ad)[1]][1]
frt_ch <- recf$input[recf$ch == min(recf$ch)[1]][1]
frt_kp <- recf$input[recf$kp == min(recf$kp)[1]][1]
frt_ks <- recf$input[recf$ks == min(recf$ks)[1]][1]
# gravels input/output rate
grt <- recg$input[recg$ks == min(recg$ks)[1]][1]
grt_ad <- recg$input[recg$ad == min(recg$ad)[1]][1]
grt_ch <- recg$input[recg$ch == min(recg$ch)[1]][1]
grt_kp <- recg$input[recg$kp == min(recg$kp)[1]][1]
grt_ks <- recg$input[recg$ks == min(recg$ks)[1]][1]

# debris-flow events remaining
der <- rec$n[rec$ks == min(rec$ks)[1]][1]
# fines events remaining
fer <- recf$n[recf$ks == min(recf$ks)[1]][1]
# gravel events remaining
ger <- recg$n[recg$ks == min(recg$ks)[1]][1]

# turnover period for debris-flow inputs
tp_ks <- 205.27607176168172
tp_ad <- 318 # ad
tp_ch <- 293 # ch
tp_kp <- 191 # kp
tp <- 208 # ks

# debris-flow events per turnover period
et <- tp * drt
et_ad <- tp_ad * drt_ad
et_ch <- tp_ch * drt_ch
et_kp <- tp_kp * drt_kp
et_ks <- tp_ks * drt_ks
# debris-flow events remaining
er_ad <- rec$n[rec$ad == min(rec$ad)[1]]
er_ch <- rec$n[rec$ch == min(rec$ch)[1]]
er_kp <- rec$n[rec$kp == min(rec$kp)[1]]
er_ks <- rec$n[rec$ks == min(rec$ks)[1]]
# er <- rec$n[rec$ch == min(rec$ch)[1]]
# er <- rec$n[rec$kp == min(rec$kp)[1]]
er <- rec$n[rec$ks == min(rec$ks)[1]]

# fluvial events per turnover period
fet <- tp * frt
fet_ad <- tp_ad * frt_ad
fet_ch <- tp_ch * frt_ch
fet_kp <- tp_kp * frt_kp
fet_ks <- tp_ks * frt_ks
# fluvial events remaining
fer_ad <- recf$n[recf$ad == min(recf$ad)[1]][1]
fer_ch <- recf$n[recf$ch == min(recf$ch)[1]][1]
fer_kp <- recf$n[recf$kp == min(recf$kp)[1]][1]
fer_ks <- recf$n[recf$ks == min(recf$ks)[1]][1]
# fer <- recf$n[recf$ch == min(recf$ch)[1]][1]
# fer <- recf$n[recf$kp == min(recf$kp)[1]][1]
fer <- recf$n[recf$ks == min(recf$ks)[1]][1]

# gravel events per turnover period
get <- tp * grt
get_ad <- tp_ad * grt_ad
get_ch <- tp_ch * grt_ch
get_kp <- tp_kp * grt_kp
get_ks <- tp_ks * grt_ks
# gravel events remaining
ger_ad <- recg$n[recg$ad == min(recg$ad)[1]][1]
ger_ch <- recg$n[recg$ch == min(recg$ch)[1]][1]
ger_kp <- recg$n[recg$kp == min(recg$kp)[1]][1]
ger_ks <- recg$n[recg$ks == min(recg$ks)[1]][1]
# ger <- recg$n[recg$ch == min(recg$ch)[1]][1]
# ger <- recg$n[recg$kp == min(recg$kp)[1]][1]
ger <- recg$n[recg$ks == min(recg$ks)[1]][1]

# total events per turnover
te <- et + fet + get
te_ad <- et_ad + fet_ad + get_ad
te_ch <- et_ch + fet_ch + get_ch
te_kp <- et_kp + fet_kp + get_kp
te_ks <- et_ks + fet_ks + get_ks

# total events remaining
tr <- er + fer + ger
tr_ad <- er_ad + fer_ad + ger_ad
tr_ch <- er_ch + fer_ch + ger_ch
tr_kp <- er_kp + fer_kp + ger_kp
tr_ks <- er_ks + fer_ks + ger_ks

# odds of flux selection
pf <- te / (te + tr)
pf_ad <- te_ad / (te_ad + tr_ad) # 0.5252573
pf_ch <- te_ch / (te_ch + tr_ch) # 0.5194105
pf_kp <- te_kp / (te_kp + tr_kp) # 0.5726348
pf_ks <- te_ks / (te_ks + tr_ks) # 0.5187506

pf_ad <- get_ad / (get_ad + ger_ad) # 0.4093586
pf_ch <- get_ch / (get_ch + ger_ch) # 0.5197172
pf_kp <- get_kp / (get_kp + ger_kp) # 0.5895043
pf_ks <- get_ks / (get_ks + ger_ks) # 0.508744
# ad = 0.7227464
# ch = 0.6812836
# kp = 0.4243889
# ks = 0.5187506

fish <- function(rt, t, k) {
  (rt*t)^k*exp(-rt*t) / factorial(k)
}


# look up optimized delivery probability and streampower coefficient k for charcoal samples
crks <- creeks_radio
crks_so <- rater1(creeks, .73, .73)
# convert lengths from to-mouth to hipchain from study area
crks_so$hip <- crks_so$ToMouth_km
crks_so$hip[crks_so$creek_name == 'knowles'] <- crks_so$hip[crks_so$creek_name == 'knowles'] - min(crks_so$hip[crks_so$creek_name == 'knowles'])

# subset bear then knowles
br <- crks_so[crks_so$creek_name == 'bear', ]
kn <- crks_so[crks_so$creek_name == 'knowles', ]

# from km to m
# switch orientation from mouth to initiation point
br$hp <- (max(br$hip) - br$hip) * 1000
kn$hp <- (max(kn$hip) - kn$hip) * 1000
# cumulative lengths
br$cm <- br$hp %>% rev %>% cumsum %>% rev
kn$cm <- kn$hp %>% rev %>% cumsum %>% rev
# upstream count
br$ct <- nrow(br):1
kn$ct <- nrow(kn):1
# average length traveled per segment
br$alt <- br$cm / br$ct
kn$alt <- kn$cm / kn$ct

# mean reach length [m]
mn_ln <- mean(c(br$hp, kn$hp))
# mean sediment distance per year [m/y]
mn_d <- mn_ln / tp

# fluvial inheritance
br$fi <- br$alt / mn_d
kn$fi <- kn$alt / mn_d


# estimate counts for fines fractions
# predict mass by width
gmod <- lm(log(wgt) ~ log(b), data = pebs[pebs$wgt > 0, ])
summary(gmod)

# combine surface pebble counts into single set
surf <- c(peb$site_1a, peb$site_1b, peb$base_1a, peb$base_1b)


# fit the in-place comminution rate
b4 <- pebs[pebs$b >= 4.5 & pebs$b <= 7.5, ]
b4c <- cdf_hook(log(b4$b))
b4d <- data.frame(b = b4$b, lb = log(b4$b), c = b4c)
b4m <- lm(c ~ lb, data = b4d)
summary(b4m)
xs <- seq(4.5, 7.5, .01)

png('comminution_fit.png',  height = 17, width = 21, units = 'cm', res = 300)
plot(b4d$b, b4d$c, log = 'x', pch = 20, col = get_palette('ocean', .1),
     xlab = 'b-axis width [mm]',
     ylab = 'CDF of sizes 4.5-7.5mm')
b4p <- b4m$coefficients[1] + (b4m$coefficients[2] * log(xs))
lines(xs, b4p, lwd = 3, col = get_palette('gold', .5))
legend('topleft', legend = c('subsurface', 'modeled'),
       fill = get_palette(c('ocean', 'gold'), .7))
text(7.3, 0.85, labels = 'R2 = 0.997')
dev.off()

# mass log-linear model from 0.25-0.6 g
gd <- pebs[pebs$wgt >= 0.25 & pebs$wgt <= .6, ]
gd$c <- cdf_hook(gd$wgt)
gdm <- lm(c ~ log(wgt), data = gd)
plot(gd$wgt, gd$c)
summary(gdm)
xm <- seq(0.25, 0.6, length.out = nrow(gd))

png('logliner_width_mass.png', height = 17, width = 21, units = 'cm', res = 300)
par(mar = c(5, 5, 5, 3))
plot(b4d$b, b4d$c, log = 'x', pch = 20, col = get_palette('ocean', .1),
     xlab = 'b-axis width [mm]',
     ylab = 'CDF of samples')
b4p <- b4m$coefficients[1] + (b4m$coefficients[2] * log(xs))
lines(xs, b4p, lwd = 2.5, col = get_palette('gold', .9))
text(7.35, 0.83, labels = 'R2 = 0.997', col = get_palette('gold', 1))
par(new = T)
plot(gd$wgt, gd$c, log = 'x', pch = 20, col = get_palette('charcoal', .1),
     xaxt = 'n', yaxt = 'n',
     xlab = '', ylab = '')
gdp <- gdm$coefficients[1] + (gdm$coefficients[2] * log(xm))
lines(xm, gdp, lwd = 2.5, col = get_palette('crimson'))
text(.58, 0.8, labels = 'R2 = 0.997', col = get_palette('crimson', 1))
axis(side = 3)
mtext('mass [g]', side = 3, line = 3)
legend('bottomright', legend = c('width', 'width fit', 'mass', 'mass fit'),
       fill = get_palette(c('ocean', 'gold', 'charcoal', 'crimson'), .7))
dev.off()

plot(emp_cdf(pebs$wgt), log ='x', xlim = c(0.15, 0.9))



# build data frames for plotting
d <- data.frame(b = surf, c = cdf_hook(surf))
ds <- data.frame(b = pebs$b, c= cdf_hook(pebs$b))

# relative rarity
d$r <- vapply(surf, function(x) rarer(x, median(surf), b4m, .05), 0)
ds$r <- vapply(ds$b, function(x) rarer(x, median(pebs$b), b4m, .05), 0)


png('rarity_fit.png', height = 17, width = 21, units = 'cm', res = 300)
plot(d$b, d$r, log ='xy', xlim = c(2, max(d$b)),
     xlab = 'b-axis width [mm]',
     ylab = 'relative rarity to D50',
     pch = 20, col = get_palette('ocean'))
points(ds$b, ds$r, pch = 20, col = get_palette('crimson'))
abline(h = 1, lty = 2, col = get_palette('charcoal', .7))
text(200, 2.5, labels = 'D50 = 1')
legend('bottomleft', legend = c('surface', 'subsurface', 'D50'),
       pch = c(20, 20, NA), lty = c(NA, NA, 2),
       col = get_palette(c('ocean', 'crimson', 'charcoal'), .77))
dev.off()

# umbrella distribution
ds$rw <- ds$b / median(pebs$b) # width relative to median
d$rw <- d$b / median(d$b) # width relative to median
# ds$rw[ds$rw < 1] <- median(pebs$b) / ds$b[ds$rw < 1]
# above median, increasing width decreases selection probability
# subtract 1 from rarity weight and max because median == 1 to get relative probability
# subtract from one to indicate decreasing selection likelihood with increasing width
ds$rw[ds$rw < 1] <- 1
ds$rw[ds$rw > 1] <- 1 - (ds$rw[ds$rw > 1] - 1) / (max(ds$rw) - 1)
d$rw[d$rw > 1] <- 1 - (d$rw[d$rw > 1] - 1) / (max(d$rw) - 1)
# below median, decreasing width decreases selection probability
# invert the width ratio
# ds$rw[ds$b < median(ds$b)] <- 1 / ds$rw[ds$b < median(ds$b)]
d$rw[d$b < median(d$b)] <- 1 / d$rw[d$b < median(d$b)]
# convert to relative probability using same transform as above median
# ds$rw[ds$rw > 1] <- 1 - (ds$rw[ds$rw > 1] - 1) / (max(ds$rw) - 1)
d$rw[d$rw > 1] <- 1 - (d$rw[d$rw > 1] - 1) / (max(d$rw) - 1)

# preferred model
m <- lm(c ~ log(b) + rw + b, data = d)
summary(m)
ms <- lm(c ~ log(b) + rw + b, data = ds)
summary(ms)

b4$c <- cdf_hook(b4$wgt)
gm <- lm(c ~ log(wgt), data = b4[b4$wgt != 0, ])
summary(gm)

dp <- predict(m, newdata = d)
msp <- predict(ms, newdata = ds)

png('gravel_fit2.png', height = 17, width = 21, units = 'cm', res = 300)
plot(sort(d$b), sort(dp), log = 'x',
     xlab = 'b-axis width [mm]', xlim = c(2, max(d$b)),
     ylab = 'CDF of samples', ylim = range(msp),
     type = 'l', lwd = 3, col = get_palette('crimson', .7))
points(emp_cdf(d$b), pch = 20, col = get_palette('charcoal', .2))
points(emp_cdf(ds$b), pch = 20, col = get_palette('hardwood', .2))
lines(sort(ds$b), sort(msp), lwd = 3, col = get_palette('crimson', .7))
legend('bottomright', legend = c('surface', 'subsurface', 'modeled'),
       fill = get_palette(c('charcoal', 'hardwood', 'crimson'), .8))
dev.off()

# estimated turnover period at selected fits

# convert pmf to cdf and calculate weighted mean transit time
dft_ad <- df_transit_ad %>% cumsum
dft_ad_mn <- weighted.mean(0:10000, df_transit_ad)
dft_ch <- df_transit_ch %>% cumsum
dft_ch_mn <- weighted.mean(0:10000, df_transit_ch)
dft_kp <- df_transit_kp %>% cumsum
dft_kp_mn <- weighted.mean(0:10000, df_transit_kp)
dft_ks <- df_transit_ks %>% cumsum
dft_ks_mn <- weighted.mean(0:10000, df_transit_ks)

# d50_mass <- exp(predict(gmod, newdata = data.frame(b = 41.15)))
d50_mass <- exp(predict(gmod, newdata = data.frame(b = median(pebs$b))))
d2_mass <- exp(predict(gmod, newdata = data.frame(b = 2)))

# sternberg mass loss coefficient over time
# estimating mass from b-axis width
# turnover period from select tests of debris-flow charcoal ages
st_ad <- -(log(d2_mass) - log(d50_mass)) / dft_ad_mn
st_ch <- -(log(d2_mass) - log(d50_mass)) / dft_ch_mn
st_kp <- -(log(d2_mass) - log(d50_mass)) / dft_kp_mn
st_ks <- -(log(d2_mass) - log(d50_mass)) / dft_ks_mn

stm <- -(log(d2_mass) - log(d50_mass)) / mn_dkn
stm <- -(log(d2_mass) - log(d50s_mass)) / mn_dkn

# sternberg equation fit to log width
stw_ad <- log(median(surf)) / dft_ad_mn
stw_ch <- log(median(surf)) / dft_ch_mn
stw_kp <- log(median(surf)) / dft_kp_mn
stw_ks <- log(median(surf)) / dft_ks_mn


xlen <- 1e4
std <- data.frame(
  ad = d50_mass * exp(-st_ad * seq(1, dft_ad_mn, length.out = xlen)),
  ch = d50_mass * exp(-st_ch * seq(1, dft_ch_mn, length.out = xlen)),
  kp = d50_mass * exp(-st_kp * seq(1, dft_kp_mn, length.out = xlen)),
  ks = d50_mass * exp(-st_ks * seq(1, dft_ks_mn, length.out = xlen)),
  st = d50_mass * exp(-stm * seq(1, dft_ad_mn, length.out = xlen)),
  oc = d50_mass * exp(-1.204 * seq(1, dft_ad_mn, length.out = xlen))
)



png('pebble_mass_loss.png', height = 17, width = 21, units = 'cm', res = 300)
magicaxis::magplot(seq(1, dft_ad_mn, length.out = xlen), std$ad, log = 'y',
     lwd = 2.5, type = 'l', col = get_palette('ocean', .7),
     xlab = 'Residence Time [y]', ylab = 'Pebble Mass [g]')
lines(seq(1, dft_ch_mn, length.out = xlen), std$ch,
      lwd = 2.5, col = get_palette('gold', .7))
lines(seq(1, dft_kp_mn, length.out = xlen), std$kp,
      lwd = 2.5, col = get_palette('violet', .7))
lines(seq(1, dft_ks_mn, length.out = xlen), std$ks,
      lwd = 2.5, col = get_palette('crimson', .7))
legend('topright', legend = c('Anderson-Darling', 'Chi-squared', 'Kuiper', 'Kolmogorov-Smirnov'),
       fill = get_palette(c('ocean', 'gold', 'violet', 'crimson'), .9))
text(50, 0.19, labels = paste0('D50s = ', round(d50_mass, 2), ' g'))
text(240, 0.012, labels = paste0('D2mm = ', round(d2_mass, 3), ' g'))
dev.off()


ss <- data.frame(g = std$ks[std$ks < 2])
ss$c <- cdf_hook(ss$g)
ssm <- lm(c ~ log(g), data = ss)
summary(ssm)

gdm
gmod
# expected median weight in grams
ew <- exp(predict(gmod, newdata = data.frame(b = median(surf))))
ews <- exp(predict(gmod, newdata = data.frame(b = median(pebs$b))))
# relative rarity
d$gr <- vapply(surf, function(x) rarer(x, ew, gdm, .05), 0)
ds$gr <- vapply(ds$b, function(x) rarer(x, ews, gdm, .05), 0) * 1e8




plot(seq(1, dft_ad_mn, length.out = 100), std$ad, log = 'y',
     lwd = 2.5, type = 'l', col = get_palette('ocean', .7),
     xlab = 'residence time [y]', ylab = 'pebble mass [g]')




std$cad <- cdf_hook(std$ad)
std$cch <- cdf_hook(std$ch)
std$cks <- cdf_hook(std$kp)
std$ckp <- cdf_hook(std$ks)

stadm <- lm(cad ~ log(ad), data = std)
stchm <- lm(cch ~ log(ch), data = std)
stksm <- lm(cks ~ log(kp), data = std)
stkpm <- lm(ckp ~ log(kp), data = std)
summary(stadm)
summary(stchm)
summary(stkpm)
summary(stksm)

log(4.5) / median(log(surf))
log(7.5) / median(log(surf))
b4m$coefficients[2] / (log(7.5) - log(4.5)) / median(log(surf))

d$rks <- vapply(surf, function(x) rarer(x, ew, stksm, .05), 0)
ds$rks <- vapply(ds$b, function(x) rarer(x, ews, stksm, .05), 0)

plot(d$b, d$gr, log = 'xy', xlim = c(2, max(d$b)),
     ylim = c(min(c(d$r, d$gr, ds$r, ds$gr)), max(c(d$r, d$gr, ds$r, ds$gr))))
points(d$b, d$r, pch = 20, col = get_palette('gold'))
points(ds$b, ds$gr, pch = 20, col = get_palette('charcoal'))
points(ds$b, ds$r, pch = 20, col = get_palette('crimson'))
plot(ds$b, ds$rks, pch = 20, col = get_palette('crimson'), log = 'xy')
plot(d$b, d$rks, pch = 20, col = get_palette('crimson'), log = 'xy')
points(d$b, d$rks, pch = 20, col = get_palette('crimson'))

df1 <- c('intercept', 'beta')
df1$coefficients <- c(1, )
gdm
off <- (median(log(surf)) * gmod$coefficients[2] + gmod$coefficients[1])
gmod
# sternberg mass loss coefficient over distance
# estimating mean travel distance from bear and knowles
# starting weight is d50 of surface gravel at knowles
stm <- log(predict(gmod, newdata = data.frame(b = median(surf)))) / (mn_ln / 1000)
# 0.816 comparable to O'Conner 2014

# mass remaining after turnover period
(mr_ad <- exp(log(sum(pebs$wgt)) - st_ad * dft_ad_mn)) # 477.58 g remain
(mr_ch <- exp(log(sum(pebs$wgt)) - st_ch * dft_ch_mn)) # 477.58 g remain
(mr_kp <- exp(log(sum(pebs$wgt)) - st_kp * dft_kp_mn)) # 477.58 g remain
(mr_ks <- exp(log(sum(pebs$wgt)) - st_ks * dft_ks_mn)) # 477.58 g remain
(mr_st <- exp(log(sum(pebs$wgt)) - stm * mn_ln / 1000)) # 477.58 g remain


# fines produced by measured gravel in turnover period [g] (50g observed)
fn <- sum(pebs$wgt) - mr_ks # 1243.835 g produced
# ratio of fines to gravel
fn / sum(pebs$wgt)
# fluvial exit probability
# % of gravel converted to fines over turnover period
fep <- fn / sum(pebs$wgt)





# bootstrap charcoal ages
s1 <- sample(dfi, length(dfi), replace = T, prob = dfp)

s2 <- s1
altered <- 0
for (i in seq_along(s1)) {
  if (runif(1) > 0.8) {
    s2[i] <- s1[i] + sample(index, 1, replace = T, prob = pex_pmf)
    altered <- altered + 1
  }
}
print(altered / length(s2))

plot(emp_cdf(s1))
plot(emp_cdf(s2))
lines(dfc)
lines(emp_cdf(charcoal$mn[charcoal$facies == 'FG']))


plot(kn, pch = '.', col = 'blue')
points(kn[kn$NODE_ID %in% crks$node_ids[crks$type == 'FG'], ])
t <- seq(0, 37, length.out = 4430)
# probability of exiting single reservoir within turnover period t (k is exit)
pex_pmf <- (fish(1-pf, t, 1)/sum(fish(1-pf,t,1)))
pex_kn <- fish(1-(pf * (mean(kn$hp) / mn_ln)), t, 1) / sum(fish(1-(pf * mean(kn$hp) / mn_ln), t, 1))
pex_br <- fish(1-(pf * (mean(br$hp) / mn_ln)), t, 1) / sum(fish(1-(pf * mean(br$hp) / mn_ln), t, 1))
pex_cd <- fish(1-(pf * (mean(crks_so$ToMouth_km[crks_so$creek_name == 'cedar']) * 1000) / mn_ln), t, 1) /
  sum(fish(1-(pf * (mean(crks_so$ToMouth_km[crks_so$creek_name == 'cedar']) * 1000) / mn_ln), t, 1))
pex_gr <- fish(1-(pf * 500 / mn_ln), t, 1) / sum(fish(1-(pf * 500 / mn_ln), t, 1))

# fluvial inheritance at select rates
# pf_ad <- 0.7227464
# pf_ch <- 0.6812836
# pf_kp <- 0.4243889
# pf_ks <- 0.5187506

tp_ad <- 318
tp_ch <- 293
tp_kp <- 191
tp_ks <- 208

sfill / tp_ad / 1416 * 0.128
sfill / tp_ch / 1515 * 0.137
sfill / tp_kp / 2864 * 0.259
sfill / tp_ks / 2109 * 0.191


t <- seq(0, 100, length.out = 4430)
pex_ad <- (fish(0.058, t, 1)/sum(fish(0.058, t, 1)))
pex_ch <- (fish(0.136, t, 1)/sum(fish(0.136, t, 1)))
pex_kp <- (fish(0.011, t, 1)/sum(fish(0.011, t, 1)))
pex_ks <- (fish(0.106, t, 1)/sum(fish(0.106, t, 1)))

fex_ad <- (fish(0.108, t, 1)/sum(fish(0.108, t, 1)))
fex_ch <- (fish(0.075, t, 1)/sum(fish(0.075, t, 1)))
fex_kp <- (fish(0.075, t, 1)/sum(fish(0.075, t, 1)))
fex_ks <- (fish(0.425, t, 1)/sum(fish(0.425, t, 1)))


png('fluvial_travel.png', height = 17, width = 21, units = 'cm', res = 300)
magicaxis::magplot(t * tp_ks , cumsum(pex_ks), log = 'x',
     xlab = 'Fluvial Travel Time [y]', xlim = c(50, 35000),
     ylab = 'CDF',
     type = 'l', lwd = 2.5, col = get_palette('crimson', .7))
lines(t*tp_ch, cumsum(pex_ch), lwd = 2.5, col = get_palette('gold', .7))
lines(t * tp_kp , cumsum(pex_kp), lwd = 2.5, col = get_palette('violet', .7))
lines(t*tp_ad, cumsum(pex_ad), lwd = 2.5, col = get_palette('ocean', .7))
lines(t * tp_ad , cumsum(fex_ad), lty = 2, lwd = 1.5, col = get_palette('ocean', .7))
lines(t * tp_ch , cumsum(fex_ch), lty = 2, lwd = 1.5, col = get_palette('gold', .7))
lines(t * tp_kp , cumsum(fex_kp), lty = 2, lwd = 1.5, col = get_palette('violet', .7))
lines(t * tp_ks , cumsum(fex_ks), lty = 2, lwd = 1.5, col = get_palette('crimson', .7))
legend('topleft', legend = c('Gravels', 'Fines', 'Anderson-Darling', 'Chi-squared', 'Kuiper', 'Kolmogorov-Smirnov'),
       col = get_palette(c('charcoal', 'charcoal', 'ocean', 'gold', 'violet', 'crimson'), .9),
       lty = c(1, 2, 1, 1, 1, 1))
dev.off()


library(magicaxis)

png('fluvial_travel.png', height = 17, width = 21, units = 'cm', res = 300)
plot(t*tp_ad, cumsum(pex_ks), log = 'x',
     xlab = 'fluvial travel time (slow fraction) [y]', xlim = c(20, 9000),
     ylab = 'CDF',
     type = 'l', lwd = 2.5, col = get_palette('ocean', .7))
lines(t*tp_ch, cumsum(pex_ks), lwd = 2.5, col = get_palette('gold', .7))
lines(t*tp_kp, cumsum(pex_ks), lwd = 2.5, col = get_palette('violet', .7))
lines(t*tp_ks, cumsum(pex_ks), lwd = 2.5, col = get_palette('crimson', .7))
legend('topleft', legend = c('Anderson-Darling', 'Chi-squared', 'Kuiper', 'Kolmogorov-Smirnov'),
       fill = get_palette(c('ocean', 'gold', 'violet', 'crimson'), .9))
dev.off()

tex <- seq(0, 1, length.out = 4430)
pex_f <- fish(fep, tex, 1)/ sum(fish(fep, tex, 1))

png('fluvial_exit.png', height = 17, width = 21, units = 'cm', res = 300)
magplot(tex * tp_ad, cumsum(pex_t), type = 'l', lwd = 2.5, col = get_palette('violet', .7),
     xlab = 'fluvial travel time (fast fraction) [y]',
     ylab = 'CDF')
lines(tex * tp_ch, cumsum(pex_t), lwd = 2.5, col = get_palette('gold', .7))
lines(tex * tp_kp, cumsum(pex_t), lwd = 2.5, col = get_palette('crimson', .7))
lines(tex * tp_ks, cumsum(pex_t), lwd = 2.5, col = get_palette('ocean', .7))
legend('topleft', legend = c('Anderson-Darling', 'Chi-squared', 'Kuiper', 'Kolmogorov-Smirnov'),
       fill = get_palette(c('ocean', 'gold', 'violet', 'crimson'), .9))
dev.off()


# gravel deposit count at study sites
# pull type from fluvial1.R
gct_kn <- nrow(crks[crks$creek_name == 'knowles' & crks$type == 'FG', ])
gct_br <- nrow(crks[crks$creek_name == 'bear' & crks$type == 'FG', ])
gct_cd <- 10
gct_gr <- 36
gct_tl <- sum(gct_kn, gct_br, gct_cd, gct_gr)

# fines deposit count at study sites
fct_kn <- nrow(crks[crks$creek_name == 'knowles' & crks$type == 'FF', ])
fct_br <- nrow(crks[crks$creek_name == 'bear' & crks$type == 'FF', ])
fct_cd <- 10
fct_gr <- 36
fct_tl <- sum(fct_kn, fct_br, fct_cd, fct_gr)

dct_kn <- nrow(crks[crks$creek_name == 'knowles' & crks$type == 'DF', ])
dct_br <- nrow(crks[crks$creek_name == 'bear' & crks$type == 'DF', ])

m_df <- (dct_br + dct_kn) * (1 - (gct_br + gct_kn) / (dct_br + dct_kn))
m_g <- fep * (gct_br + gct_kn)
# fines exit probability
ffep <- (m_df + m_g) / (m_df + m_g + fct_br + fct_kn)
# remainder probability (fines)
rpf <- 1 - ffep

# fines exit pmfs
ftrap <- 4
fex_t <- fish(ffep^ftrap, t, 1) / sum(fish(ffep^ftrap, t, 1))
# fpex_t <- fish(ffep, tex, 1)/ sum(fish(ffep, tex, 1))

# debris-flow pmf and index
index <- char_pmfs %>% rownames %>% as.numeric %>% rev
dfc <- emp_cdf(charcoal$mn[charcoal$facies == 'DF'])
dfi <- dfc[ , 1]
dfp <- to_pmf(dfc[ , 2])



#
dfi <- dfi + 50
hill_source <- dfi[dfi < tp_ks]
hill_pmf <- dfp[dfi < tp_ks] / sum(dfp[dfi < tp_ks])
plot(hill_source, cumsum(hill_pmf))

gr_fit <- gravel_fit(n = 10, win = .33)
gr_ad <- gr_fit[[1]]
gr_ch <- gr_fit[[2]]
gr_kp <- gr_fit[[3]]
gr_ks <- gr_fit[[4]]

fn_fit <- fine_fit(mod = .4, hill_rt = 0.0)
fn_ad <- fn_fit[[1]]
fn_ch <- fn_fit[[2]]
fn_kp <- fn_fit[[3]]
fn_ks <- fn_fit[[4]]

png('gravel_fits.png', height = 17, width = 21, units = 'cm', res = 300)
magplot(emp_cdf(charcoal$mn[charcoal$facies == 'FF']), pch = 20, col = get_palette('coral'),
        xlab = 'Charcoal Age [y]', ylab = 'CDF', xlim = c(0, 17000))
points(dfc, pch = 20, col = get_palette('hardwood'))
points(emp_cdf(charcoal$mn[charcoal$facies == 'FG']), pch = 20, col = get_palette('charcoal'))
lines(gr_ad, lwd = 2.5, col = get_palette('ocean'))
lines(gr_ch, lwd = 2.5, col = get_palette('gold'))
lines(gr_kp, lwd = 2.5, col = get_palette('violet'))
lines(gr_ks, lwd = 2.5, col = get_palette('crimson'))
# lines(gf[[1]], lwd = 2.5, col = get_palette('ocean'))
# lines(gf[[2]], lwd = 2.5, col = get_palette('gold'))
# lines(gf[[3]], lwd = 2.5, col = get_palette('violet'))
# lines(gf[[4]], lwd = 2.5, col = get_palette('crimson'))
legend('bottomright', legend = c('Anderson-Darling', 'Chi-squared', 'Kuiper', 'Kolmogorov-Smirnov', 'Debris Flows', 'Gravels', 'Fines'),
       pch = c(rep(NA, 4), rep(20, 3)), lty = c(rep(1, 4), rep(NA, 3)),
       col = get_palette(c('ocean', 'gold', 'violet', 'crimson', 'hardwood', 'charcoal', 'coral'), .77))
dev.off()

png('fines_fits.png', height = 17, width = 21, units = 'cm', res = 300)
magplot(emp_cdf(charcoal$mn[charcoal$facies == 'FF']), pch = 20, col = get_palette('coral'),
        xlab = 'charcoal age', ylab = 'CDF of samples', xlim = c(0, 17000))
points(dfc, pch = 20, col = get_palette('hardwood'))
points(emp_cdf(charcoal$mn[charcoal$facies == 'FG']), pch = 20, col = get_palette('charcoal'))
lines(fn_ad, lwd = 2.5, col = get_palette('ocean'))
lines(fn_ch, lwd = 2.5, col = get_palette('gold'))
lines(fn_kp, lwd = 2.5, col = get_palette('violet'))
lines(fn_ks, lwd = 2.5, col = get_palette('crimson'))
legend('bottomright', legend = c('anderson-darling', 'chi-squared', 'kuiper', 'kolmogorov-smirnov', 'debris flows', 'gravels', 'fines'),
       pch = c(rep(NA, 4), rep(20, 3)), lty = c(rep(1, 4), rep(NA, 3)),
       col = get_palette(c('ocean', 'gold', 'violet', 'crimson', 'hardwood', 'charcoal', 'coral'), .77))
dev.off()

df_test <- data.table::fread('/home/erik/output/df_test.csv')
df_test <- unlist(df_test)


lines(emp_cdf(gr_test), lwd = 2)

lines(emp_cdf(df_test), lwd = 2)

library(plot3D)
gr_test <- data.table::fread('/home/erik/output/gravels_ks_1kx_1000.csv')
min(gr_test$ks2)
(gr_ks_flux <- gr_test$capture_rate_gravels[gr_test$ks2 == min(gr_test$ks2)][1])
(gr_ks_stor <- gr_test$storage_rate_gravels[gr_test$ks2 == min(gr_test$ks2)][1])
plot3D::points3D(gr_test$capture_rate_gravels, gr_test$storage_rate_gravels, gr_test$ks2, ticktype = 'detailed',
                 xlab = 'flux rate', ylab = 'storage rate', zlab = 'ks test', pch = 20,
                 phi = 30, theta = 310)
plot(emp_cdf(gr_test$ks2))
magplot(gr_test$storage_rate, gr_test$ks2,
        xlab = 'storage rate', ylab = 'ks test')
plot3D::points3D(gr_test$capture_rate_fines, gr_test$storage_rate_fines, gr_test$ks2, ticktype = 'detailed',
                 xlab = 'flux rate', ylab = 'storage rate', zlab = 'ks test', pch = 20,
                 phi = 30, theta = 310)

gr_test_ad <- data.table::fread('/home/erik/output/gravels_ad_1kx_1000.csv')
min(gr_test_ad$ad2)
(gr_ad_flux <- gr_test_ad$capture_rate_gravels[gr_test_ad$ad2 == min(gr_test_ad$ad2)][1])
(gr_ad_stor <- gr_test_ad$storage_rate_gravels[gr_test_ad$ad2 == min(gr_test_ad$ad2)][1])
plot3D::points3D(gr_test_ad$capture_rate_gravels, gr_test_ad$storage_rate_gravels, gr_test_ad$ad2, ticktype = 'detailed',
                 xlab = 'capture rate', ylab = 'storage rate', zlab = 'ad test', pch = 20,
                 phi = 20, theta = 310)
gr_test_ad <- gr_test_ad[gr_test_ad$ad2 < 2, ]
plot3D::points3D(gr_test_ad$capture_rate_gravels, gr_test_ad$storage_rate_gravels, gr_test_ad$ad2, ticktype = 'detailed',
                 xlab = 'capture rate', ylab = 'storage rate', zlab = 'ad test', pch = 20,
                 phi = 20, theta = 310)

plot(emp_cdf(gr_test_ad$ad2))
(gr_ad_flux <- gr_test_ad$capture_rate_fines[gr_test_ad$ad2 == min(gr_test_ad$ad2)][1])
(gr_ad_stor <- gr_test_ad$storage_rate_fines[gr_test_ad$ad2 == min(gr_test_ad$ad2)][1])
plot3D::points3D(gr_test_ad$capture_rate_fines, gr_test_ad$storage_rate_fines, gr_test_ad$ad2, ticktype = 'detailed',
                 xlab = 'capture rate', ylab = 'storage rate', zlab = 'ad test', pch = 20)
magplot(gr_test_ad$storage_rate, gr_test_ad$ad2,
        xlab = 'storage rate', ylab = 'Anderson-Darling test',
        pch = 20, col = get_palette('crimson'))

gr_test_ch <- data.table::fread('/home/erik/output/gravels_ch_1kx_1000.csv')
gr_test_ch <- rbind(gr_test_ch,
                    data.table::fread('/home/erik/output/gravels_ch_10kx_1001.csv'))
plot(emp_cdf(gr_test_ch$ch))
min(gr_test_ch$ch)
(gr_ch_flux <- gr_test_ch$capture_rate_gravels[gr_test_ch$ch == min(gr_test_ch$ch)][1])
(gr_ch_stor <- gr_test_ch$storage_rate_gravels[gr_test_ch$ch == min(gr_test_ch$ch)][1])
plot3D::points3D(gr_test_ch$capture_rate_gravels, gr_test_ch$storage_rate_gravels, gr_test_ch$ch, ticktype = 'detailed',
                 xlab = 'capture rate', ylab = 'storage rate', zlab = 'ch test', pch = 20,
                 phi = 30, theta = 310)
gr_test_ch <- gr_test_ch[gr_test_ch$ch < 4]
plot3D::points3D(gr_test_ch$capture_rate_gravels, gr_test_ch$storage_rate_gravels, gr_test_ch$ch, ticktype = 'detailed',
                 xlab = 'capture rate', ylab = 'storage rate', zlab = 'ch test', pch = 20,
                 phi = 30, theta = 310)

par(mar = c(5, 5, 3, 5))
magplot(gr_test_ch$storage_rate, gr_test_ch$ch,
        xlab = 'storage rate', ylab = 'Chi-squared test',
        pch = 20, col = get_palette('gold'))
par(new = T)
magplot(gr_test_ad$storage_rate, gr_test_ad$ad2,
        pch = 20, col = get_palette('crimson'),
        xaxt = 'n', yaxt = 'n',
        xlab = '', ylab = '')
axis(side = 4)
mtext('Anderson-Darling test', side = 4, line = 3)

gr_test_kp <- data.table::fread('/home/erik/output/gravels_kp_1kx_1000.csv')
min(gr_test_kp$kp)
(gr_kp_flux <- gr_test_kp$capture_rate_gravels[gr_test_kp$kp == min(gr_test_kp$kp)][1])
(gr_kp_stor <- gr_test_kp$storage_rate_gravels[gr_test_kp$kp == min(gr_test_kp$kp)][1])
plot3D::points3D(gr_test_kp$capture_rate_gravels, gr_test_kp$storage_rate_gravels, gr_test_kp$kp, ticktype = 'detailed',
                 xlab = 'capture rate', ylab = 'storage rate', zlab = 'kp test', pch = 20,
                 phi = 30, theta = 310)
gr_test_kp <- gr_test_kp[gr_test_kp$kp < .3, ]
plot3D::points3D(gr_test_kp$capture_rate_gravels, gr_test_kp$storage_rate_gravels, gr_test_kp$kp, ticktype = 'detailed',
                 xlab = 'capture rate', ylab = 'storage rate', zlab = 'kp test', pch = 20,
                 phi = 30, theta = 310)

(gr_kp_fluxf <- gr_test_kp$capture_rate_fines[gr_test_kp$kp == min(gr_test_kp$kp)][1])
(gr_kp_storf <- gr_test_kp$storage_rate_fines[gr_test_kp$kp == min(gr_test_kp$kp)][1])
plot3D::points3D(gr_test_kp$capture_rate_fines, gr_test_kp$storage_rate_fines, gr_test_kp$kp, ticktype = 'detailed',
                 xlab = 'capture rate', ylab = 'storage rate', zlab = 'kp test', pch = 20)
magplot(gr_test_kp$storage_rate, gr_test_kp$kp,
        xlab = 'storage rate', ylab = 'Kuiper test',
        pch = 20, col = get_palette('violet'))


log(gr_ks_stor) / log(1 - gr_ks_flux)
log(gr_ks_stor) / log(1 - pf_ks)
log(gr_ad_stor) / log(1 - gr_ad_flux)
log(gr_ad_stor) / log(1 - pf_ad)
log(gr_ch_stor) / log(1 - gr_ch_flux)
log(gr_ch_stor) / log(1 - pf_ch)
log(gr_kp_stor) / log(gr_kp_flux)
log(gr_kp_stor) / log(0.48)

library(magicaxis)

gr_stereo <- data.table::fread('/home/erik/output/gravels_stereo.csv')
gr_stereo_ad <- data.table::fread('/home/erik/output/gravels_stereo_ad.csv')
gr_stereo_ch <- data.table::fread('/home/erik/output/gravels_stereo_ch.csv')
gr_stereo_kp <- data.table::fread('/home/erik/output/gravels_stereo_kp.csv')
gr_stereo_ks <- data.table::fread('/home/erik/output/gravels_stereo_ks.csv')
fn_stereo_ad <- data.table::fread('/home/erik/output/fines_stereo_ad.csv')
fn_stereo_kp <- data.table::fread('/home/erik/output/fines_stereo_kp.csv')
magicaxis::magplot(emp_cdf(charcoal$mn[charcoal$facies == 'FF']), pch = 20, col = get_palette('coral'),
        xlab = 'charcoal age', ylab = 'CDF of samples', xlim = c(0, 17000))
points(dfc, pch = 20, col = get_palette('hardwood'))
points(emp_cdf(charcoal$mn[charcoal$facies == 'FG']), pch = 20, col = get_palette('charcoal'))
lines(emp_cdf(fn_stereo_kp %>% unlist), lwd = 1.5, col = get_palette('violet', .7))
lines(emp_cdf(fn_stereo_ad %>% unlist), lwd = 1.5, col = get_palette('ocean', .7))

lines(emp_cdf(gr_stereo_ks %>% unlist), lwd = 1.5, col = get_palette('crimson', .7))
lines(emp_cdf(gr_stereo_ad %>% unlist), lwd = 1.5, col = get_palette('ocean', .7))
lines(emp_cdf(gr_stereo_ch %>% unlist), lwd = 1.5, col = get_palette('gold', .7))
lines(emp_cdf(gr_stereo_kp %>% unlist), lwd = 1.5, col = get_palette('violet', .7))

df_stereo <- data.table::fread('/home/erik/output/debris_flows_stereo.csv')
lines(emp_cdf(df_stereo %>% unlist))







# bootstrap gravel ages 3
gravel_synth <- function(win = fep) {
  gad <- sample(dfi, length(dfi), replace = TRUE, prob = dfp)
  gch <- sample(dfi, length(dfi), replace = TRUE, prob = dfp)
  gkp <- sample(dfi, length(dfi), replace = TRUE, prob = dfp)
  gks <- sample(dfi, length(dfi), replace = TRUE, prob = dfp)

  for (i in 1:length(dfi)) {
    rd_ad <- runif(1)
    rd_ch <- runif(1)
    rd_kp <- runif(1)
    rd_ks <- runif(1)

    # pf is exit probability
    # rd > pf means deposit did not exit
    # pex controls steepness of curve along x-axis
    # win controls % of CDF affected on y-axis
    if (rd_ad <= (1 - pf_ad) * win) {
      gad[i] <- gad[i] + sample(index, 1, prob = pex_ad)
    }
    if (rd_ch <= (1 - pf_ch) * win) {
      gch[i] <- gch[i] + sample(index, 1, replace = T, prob = pex_ch)
    }
    if (rd_kp <= (1 - pf_kp) * win) {
      gkp[i] <- gkp[i] + sample(index, 1, replace = T, prob = pex_kp)
    }
    if (rd_ks <= (1 - pf_ks) * win) {
      gks[i] <- gks[i] + sample(index, 1, replace = T, prob = pex_ks)
    }

  }
  list(gad, gch, gkp, gks)
}
# gravel fit for gravel synth 3
gravel_fit <- function(n = 10, win = fep) {
  fits <- gravel_synth(win)
  gad <- fits[[1]]
  gch <- fits[[2]]
  gkp <- fits[[3]]
  gks <- fits[[4]]
  for (i in 2:n) {
    fits <- gravel_synth(win)
    gad <- c(gad, fits[[1]])
    gch <- c(gch, fits[[2]])
    gkp <- c(gkp, fits[[3]])
    gks <- c(gks, fits[[4]])
  }
  list(emp_cdf(gad), emp_cdf(gch), emp_cdf(gkp), emp_cdf(gks))
}


# bootstrap fines ages 3
fines_synth <- function(mod = ffep, win = fep, n = 10) {
  df <- charcoal$mn[charcoal$facies == 'DF'] # debris flow age obs
  dfc <- emp_cdf(df) # cdf of ages
  dfi <- dfc[ , 1] # index of cdf ages
  dfp <- to_pmf(dfc[ , 2]) # pmf of index ages
  gf <- gravel_fit(n = n, win = win)
  fad <- 0 # synth deposits
  fch <- 0
  fkp <- 0
  fks <- 0
  for (i in 1:length(dfi)) {
    fad[i] <- sample(c(sample(dfi, 135, replace = T, prob = dfp), sample(gf[[1]][,1], 83, replace = T, prob = to_pmf(gf[[1]][,2]))), 1)
    fch[i] <- sample(c(sample(dfi, 135, replace = T, prob = dfp), sample(gf[[2]][,1], 83, replace = T, prob = to_pmf(gf[[2]][,2]))), 1)
    fkp[i] <- sample(c(sample(dfi, 135, replace = T, prob = dfp), sample(gf[[3]][,1], 83, replace = T, prob = to_pmf(gf[[3]][,2]))), 1)
    fks[i] <- sample(c(sample(dfi, 135, replace = T, prob = dfp), sample(gf[[4]][,1], 83, replace = T, prob = to_pmf(gf[[4]][,2]))), 1)
    # fad[i] <- sample(c(rep(sample(dfi, 1, prob = dfp), 3), sample(gf[[1]], 1)), 1)
    # fch[i] <- sample(c(rep(sample(dfi, 1, prob = dfp), 3), sample(gf[[2]], 1)), 1)
    # fkp[i] <- sample(c(rep(sample(dfi, 1, prob = dfp), 3), sample(gf[[3]], 1)), 1)
    # fks[i] <- sample(c(rep(sample(dfi, 1, prob = dfp), 3), sample(gf[[4]], 1)), 1)

    rd_ad <- runif(1)
    rd_ch <- runif(1)
    rd_kp <- runif(1)
    rd_ks <- runif(1)

    # pf is exit probability
    # rd > pf means deposit did not exit
    if (rd_ad <= (1 - ffep) * mod) {
      fad[i] <- fad[i] + sample(index, 1, prob = fex_t)
    }
    if (rd_ch <= (1 - ffep) * mod) {
      fch[i] <- fch[i] + sample(index, 1, prob = fex_t)
    }
    if (rd_kp <= (1 - ffep) * mod) {
      fkp[i] <- fkp[i] + sample(index, 1, prob = fex_t)
    }
    if (rd_ks <= (1 - ffep) * mod) {
      fks[i] <- fks[i] + sample(index, 1, prob = fex_t)
    }

  }
  list(fad, fch, fkp, fks)
}
# fines fit for fines synth 3
fine_fit <- function(mod = ffep, win = fep, n = 10, hill_rt = 0.1) {
  fits <- fines_synth(mod, win, n)
  gad <- fits[[1]]
  gch <- fits[[2]]
  gkp <- fits[[3]]
  gks <- fits[[4]]
  for (i in 2:n) {
    fits <- fines_synth(mod, win, n)
    gad <- c(gad, fits[[1]])
    gch <- c(gch, fits[[2]])
    gkp <- c(gkp, fits[[3]])
    gks <- c(gks, fits[[4]])
  }
  list(emp_cdf(c(gad, hills(gad, hill_rt))),
       emp_cdf(c(gch, hills(gch, hill_rt))),
       emp_cdf(c(gkp, hills(gkp, hill_rt))),
       emp_cdf(c(gks, hills(gks, hill_rt))))
}

hills <- function(obs, rate = 0.1) {
  n <- round(length(obs) * rate)
  hills <- sample(hill_source, n, replace = T, prob = hill_pmf)
  c(hill_source, hills)
}

hills()




png('fines_travel.png', height = 17, width = 21, units = 'cm', res = 300)
plot(t*tp_ad, cumsum(fex_t), log = 'x',
     xlab = 'fines travel time (slow fraction) [y]', xlim = c(tp_kp, 12000),
     ylab = 'CDF',
     type = 'l', lwd = 2.5, col = get_palette('ocean', .7))
lines(t*tp_ch, cumsum(fex_t), lwd = 2.5, col = get_palette('gold', .7))
lines(t*tp_kp, cumsum(fex_t), lwd = 2.5, col = get_palette('violet', .7))
lines(t*tp_ks, cumsum(fex_t), lwd = 2.5, col = get_palette('crimson', .7))
legend('topleft', legend = c('Anderson-Darling', 'Chi-squared', 'Kuiper', 'Kolmogorov-Smirnov'),
       fill = get_palette(c('ocean', 'gold', 'violet', 'crimson'), .9))
dev.off()

png('fines_exit.png', height = 17, width = 21, units = 'cm', res = 300)
magplot(tex * tp_ad, cumsum(fpex_t), type = 'l', lwd = 2.5, col = get_palette('violet', .7),
     xlab = 'fines travel time (fast fraction) [y]',
     ylab = 'CDF')
lines(tex * tp_ch, cumsum(fpex_t), lwd = 2.5, col = get_palette('gold', .7))
lines(tex * tp_kp, cumsum(fpex_t), lwd = 2.5, col = get_palette('crimson', .7))
lines(tex * tp_ks, cumsum(fpex_t), lwd = 2.5, col = get_palette('ocean', .7))
legend('topleft', legend = c('Anderson-Darling', 'Chi-squared', 'Kuiper', 'Kolmogorov-Smirnov'),
       fill = get_palette(c('ocean', 'gold', 'violet', 'crimson'), .9))
dev.off()


fex_kn <- fish(1-(fep * (mean(kn$hp) / mn_ln)), t, 1) / sum(fish(1-(fep * mean(kn$hp) / mn_ln), t, 1))
fex_br <- fish(1-(fep * (mean(br$hp) / mn_ln)), t, 1) / sum(fish(1-(fep * mean(br$hp) / mn_ln), t, 1))
fex_cd <- fish(1-(fep * (mean(crks_so$ToMouth_km[crks_so$creek_name == 'cedar']) * 1000) / mn_ln), t, 1) /
  sum(fish(1-(fep * (mean(crks_so$ToMouth_km[crks_so$creek_name == 'cedar']) * 1000) / mn_ln), t, 1))
fex_gr <- fish(1-(fep * 500 / mn_ln), t, 1) / sum(fish(1-(fep * 500 / mn_ln), t, 1))


png('fluvial_traversal.png', height = 17, width = 21, units = 'cm', res = 300)
magplot(t*tp_ad, cumsum(pex_ad), log = 'x',
     xlab = 'Fluvial Traversal Time [y]', xlim = c(10, 9000),
     ylab = 'CDF',
     type = 'l', lwd = 2.5, col = get_palette('ocean', .7))
lines(t*tp_ch, cumsum(pex_ch), lwd = 2.5, col = get_palette('gold', .7))
lines(t*tp_kp, cumsum(pex_kp), lwd = 2.5, col = get_palette('violet', .7))
lines(t*tp_ks, cumsum(pex_ks), lwd = 2.5, col = get_palette('crimson', .7))
lines(t*tp_ad, cumsum(fex_t), lwd = 1.75, lty = 2, col = get_palette('ocean', .7))
lines(t*tp_ch, cumsum(fex_t), lwd = 1.75, lty = 2, col = get_palette('gold', .7))
lines(t*tp_kp, cumsum(fex_t), lwd = 1.75, lty = 2, col = get_palette('violet', .7))
lines(t*tp_ks, cumsum(fex_t), lwd = 1.75, lty = 2, col = get_palette('crimson', .7))
legend('topleft', legend = c('Anderson-Darling', 'Chi-squared', 'Kuiper', 'Kolmogorov-Smirnov',
                             'Gravels', 'Fines'),
       col = get_palette(c('ocean', 'gold', 'violet', 'crimson', 'charcoal', 'charcoal'), .9),
       lty = c(rep(NA, 4), 1, 2), pch = c(rep(20, 4), NA, NA))
dev.off()

png('fast_fraction.png', height = 17, width = 21, units = 'cm', res = 300)
magplot(tex * tp_ad, cumsum(pex_f), type = 'l', lwd = 2.5, col = get_palette('ocean', .7),
        xlab = 'fluvial travel time (fast fraction) [y]',
        ylab = 'CDF')
lines(tex * tp_ch, cumsum(pex_f), lwd = 2.5, col = get_palette('gold', .7))
lines(tex * tp_kp, cumsum(pex_f), lwd = 2.5, col = get_palette('violet', .7))
lines(tex * tp_ks, cumsum(pex_f), lwd = 2.5, col = get_palette('crimson', .7))
lines(tex * tp_ad, cumsum(fpex_t), lwd = 1.75, lty = 2, col = get_palette('ocean', .7))
lines(tex * tp_ch, cumsum(fpex_t), lwd = 1.75, lty = 2, col = get_palette('gold', .7))
lines(tex * tp_kp, cumsum(fpex_t), lwd = 1.75, lty = 2, col = get_palette('violet', .7))
lines(tex * tp_ks, cumsum(fpex_t), lwd = 1.75, lty = 2, col = get_palette('crimson', .7))
legend('topleft', legend = c('Anderson-Darling', 'Chi-squared', 'Kuiper', 'Kolmogorov-Smirnov',
                             'Gravels', 'Fines'),
       col = get_palette(c('ocean', 'gold', 'violet', 'crimson', 'charcoal', 'charcoal'), .9),
       lty = c(rep(NA, 4), 1, 2), pch = c(rep(20, 4), NA, NA))
dev.off()


# convert pmf to cdf and calculate weighted mean transit time
dft_ad <- df_transit_ad %>% cumsum
dft_ad_mn <- weighted.mean(0:10000, df_transit_ad)
dft_ch <- df_transit_ch %>% cumsum
dft_ch_mn <- weighted.mean(0:10000, df_transit_ch)
dft_kp <- df_transit_kp %>% cumsum
dft_kp_mn <- weighted.mean(0:10000, df_transit_kp)
dft_ks <- df_transit_ks %>% cumsum
dft_ks_mn <- weighted.mean(0:10000, df_transit_ks)

#medians
sum((dft_ad > 0.5 & dft_ad < 0.501) * c(0:10000))
sum((dft_ch > 0.5 & dft_ch < 0.501) * c(0:10000))
sum((dft_kp > 0.5 & dft_kp < 0.501) * c(0:10000))
sum((dft_ks > 0.5 & dft_ks < 0.501) * c(0:10000))

png('df_transit_times.png', height = 17, width = 21, units = 'cm', res = 300)
magplot(dft_ad, type = 'l', lwd = 2, col = get_palette('ocean'), xlim = c(1, 3000),
        xlab = 'transit time [y]', ylab = 'CDF of deposits', log = 'x')
abline(v = dft_ad_mn, lwd = 2, lty = 2, col = get_palette('ocean'))
lines(dft_ch, type = 'l', lwd = 2, col = get_palette('gold'))
abline(v = dft_ch_mn, lwd = 2, lty = 2, col = get_palette('gold'))
lines(dft_kp, type = 'l', lwd = 2, col = get_palette('violet'))
abline(v = dft_kp_mn, lwd = 2, lty = 2, col = get_palette('violet'))
lines(dft_ks, type = 'l', lwd = 2, col = get_palette('crimson'))
abline(v = dft_ks_mn, lwd = 2, lty = 2, col = get_palette('crimson'))
legend('topleft', legend = c('anderson-darling', 'chi-squared', 'kuiper',
                             'kolmogorov-smirnov', 'mean'),
       col = get_palette(c('ocean', 'gold', 'violet', 'crimson', 'charcoal'), .8),
       lty = c(rep(1, 4), 2), lwd = 2)
text(500, 0.07, labels = '191', col = get_palette('violet', .7))
text(500, 0.14, labels = '208', col = get_palette('crimson', .7))
text(500, 0.21, labels = '293', col = get_palette('gold', .7))
text(500, 0.28, labels = '318', col = get_palette('ocean', .7))
dev.off()



plot(index, cumsum(fex_kn), log = 'x')
lines(index, cumsum(fex_br))
lines(index, cumsum(fex_cd))
lines(index, cumsum(fex_gr))

sample(c(gct_kn, gct_br, gct_cd, gct_gr), 10, replace = T, prob = c(gct_kn, gct_br, gct_cd, gct_gr) / gct_tl)
sample(1:4, 20, replace = T, prob = c(gct_kn, gct_br, gct_cd, gct_gr) / gct_tl)


pr <- sample(list(pex_kn, pex_br, pex_cd, pex_gr), 20, replace = T, prob = c(gct_kn, gct_br, gct_cd, gct_gr) / gct_tl)[[1]]


# bootstrap gravel ages
s1 <- sample(dfi, length(dfi), replace = T, prob = dfp)

s2 <- s1
altered <- 0
for (i in seq_along(s1)) {
  if (runif(1) > .8) {
    pr <- sample(list(pex_kn, pex_br, pex_cd, pex_gr), 1, replace = T, prob = c(gct_kn, gct_br, gct_cd, gct_gr) / gct_tl)[[1]]
    s2[i] <- s1[i] + sample(index, 1, replace = T, prob = pr)
    altered <- altered + 1
  }
}
print(altered / length(s2))

plot(emp_cdf(s2))
lines(dfc)
lines(emp_cdf(charcoal$mn[charcoal$facies == 'FG']))

# bootstrap fines ages

s3 <- list(s1, s2)
altered <- 0
s4 <- 0
for (i in 1:length(s1)) {
  s4[i] <- sample(sample(s3)[[1]], 1)
  if (runif(1) > .8) {
    pr <- sample(list(fex_kn, fex_br, fex_cd, fex_gr), 1, replace = T, prob = c(gct_kn, gct_br, gct_cd, gct_gr) / gct_tl)[[1]]
    s4[i] <- s4[i] + sample(index, 1, replace = T, prob = pr)
    altered <- altered + 1
  }
}
print(altered / length(s4))

plot(emp_cdf(s4), type = 'l', lwd = 3, col = get_palette('gold'))
points(dfc, pch = 20, col = get_palette('hardwood'))
points(emp_cdf(charcoal$mn[charcoal$facies == 'FF']), pch = 20, col = get_palette('gold'))
points(emp_cdf(charcoal$mn[charcoal$facies == 'FG']), pch = 20, col = get_palette('ocean'))
lines(emp_cdf(s2), col = get_palette('ocean'), lwd = 3)

crks_so$creek_name %>% factor %>% levels

index <- char_pmfs %>% rownames %>% as.numeric %>% rev
plot(index, cumsum(pex_kn), log = 'x')
lines(index, cumsum(pex_br))
lines(index, cumsum(pex_cd))
lines(index, cumsum(pex_gr))

cfl <- convo_add(pex_pmf, df_pmf, index)
tfl <- df_pmf * (te/(te+tr)) + cfl * (tr/(te+tr))
png('predicted_fluvial_age.png', height = 20, width = 24, units = 'cm', res = 300)
plot(index, tfl %>% to_cdf, type = 'l', col = get_palette('crimson'), lwd = 3,
     xlab = 'charcoal age [years]',
     ylab = 'CDF', log = 'x', ylim = c(0,1))
lines(emp_cdf(charcoal$mn[charcoal$facies == 'FF']), col = get_palette('sky'), lwd = 3)
lines(emp_cdf(charcoal$mn[charcoal$facies == 'FG']), col = get_palette('ocean'), lwd = 3)
lines(emp_cdf(charcoal$mn[charcoal$facies == 'DF']), col = get_palette('hardwood'), lwd = 3)
lines(emp_cdf(df_stereo), col = get_palette('crimson'), lwd = 3)
legend('bottomright', legend = c('gravels', 'fines', 'predicted'),
       fill = get_palette(c('ocean', 'sky', 'crimson'), .8))
dev.off()

# probability of exit before turnover period
tb <- seq(0, 205, length.out = 4430)
pex_b <- (fish(pf/tp, tb, 1)/sum(fish(pf/tp,tb,1)))
plot(tb , pex_b)
plot(tb , cumsum(pex_b))
plot(t * tp, pex_pmf)
lines(tb * tp, pex_b)
sum(pex_b)



bpex <- sum(pex_b[tb < dfi[1]])
for (i in 2:length(dfi)) {
  bpex[i] <- sum(pex_b[tb > dfi[i-1] & pex_b <= dfi[i]])
}
bpex <- bpex / sum(bpex)
plot(dfi, bpex, log = 'x')
sample(dfi, 10, prob = bpex)

# probability of exit before turnover period
fex_b <- (fish(fep/tp, tb, 1)/sum(fish(fep/tp,tb,1)))
bfex <- sum(fex_b[tb < dfi[1]])
for (i in 2:length(dfi)) {
  bfex[i] <- sum(fex_b[tb > dfi[i-1] & fex_b <= dfi[i]])
}
bfex <- bfex / sum(bfex)
plot(dfi, bfex, log = 'x')
sample(dfi, 10, prob = bfex)

# fluvial inheritance at select rates
pf_ad <- 0.7227464
pf_ch <- 0.6812836
pf_kp <- 0.4243889
pf_ks <- 0.5187506

# ad = 0.7227464
# ch = 0.6812836
# kp = 0.4243889
# ks = 0.5187506


# fit select flux rates to study area contr area
# fluvial gravel observations at each study area (for selection probability)
gct_kn <- nrow(crks[crks$creek_name == 'knowles' & crks$type == 'FG', ])
gct_br <- nrow(crks[crks$creek_name == 'bear' & crks$type == 'FG', ])
gct_cd <- 10
gct_gr <- 36
gct_tl <- sum(gct_kn, gct_br, gct_cd, gct_gr)
gct_p <- c(gct_kn, gct_br, gct_cd, gct_gr) / gct_tl

fct_kn <- nrow(crks[crks$creek_name == 'knowles' & crks$type == 'FF', ])
fct_br <- nrow(crks[crks$creek_name == 'bear' & crks$type == 'FF', ])

# contributing area weight for selection probability
caw_kn <- mean(kn$hp) / mn_ln
caw_br <- mean(br$hp) / mn_ln
caw_cd <- mean(crks_so$ToMouth_km[crks_so$creek_name == 'cedar']) * 1000 / mn_ln
caw_gr <- 500 / mn_ln

# anderson-darling
# fines exit pmfs (slow fraction)
fex_ad_kn <- fish(1-(pf_ad * caw_kn), t, 1) / sum(fish(1-(pf_ad * caw_kn), t, 1))
fex_ad_br <- fish(1-(pf_ad * caw_br), t, 1) / sum(fish(1-(pf_ad * caw_br), t, 1))
fex_ad_cd <- fish(1-(pf_ad * caw_cd), t, 1) / sum(fish(1-(pf_ad * caw_cd), t, 1))
fex_ad_gr <- fish(1-(pf_ad * caw_gr), t, 1) / sum(fish(1-(pf_ad * caw_gr), t, 1))
# gravels (slow fraction)
gex_ad_kn <- fish((pf_ad * caw_kn), t, 1) / sum(fish((pf_ad * caw_kn), t, 1))
gex_ad_br <- fish((pf_ad * caw_br), t, 1) / sum(fish((pf_ad * caw_br), t, 1))
gex_ad_cd <- fish((pf_ad * caw_cd), t, 1) / sum(fish((pf_ad * caw_cd), t, 1))
gex_ad_gr <- fish((pf_ad * caw_gr), t, 1) / sum(fish((pf_ad * caw_gr), t, 1))

fex_ad_p <- list(fex_ad_kn, fex_ad_br, fex_ad_cd, fex_ad_gr)
gex_ad_p <- list(gex_ad_kn, gex_ad_br, gex_ad_cd, gex_ad_gr)

# chi-squared
# fines exit pmfs (slow fraction)
fex_ch_kn <- fish(1-(pf_ad * caw_kn), t, 1) / sum(fish(1-(pf_ad * caw_kn), t, 1))
fex_ch_br <- fish(1-(pf_ad * caw_br), t, 1) / sum(fish(1-(pf_ad * caw_br), t, 1))
fex_ch_cd <- fish(1-(pf_ad * caw_cd), t, 1) / sum(fish(1-(pf_ad * caw_cd), t, 1))
fex_ch_gr <- fish(1-(pf_ad * caw_gr), t, 1) / sum(fish(1-(pf_ad * caw_gr), t, 1))

# gravels (slow fraction)
gex_ch_kn <- fish((pf_ch * caw_kn), t, 1) / sum(fish((pf_ch * caw_kn), t, 1))
gex_ch_br <- fish((pf_ch * caw_br), t, 1) / sum(fish((pf_ch * caw_br), t, 1))
gex_ch_cd <- fish((pf_ch * caw_cd), t, 1) / sum(fish((pf_ch * caw_cd), t, 1))
gex_ch_gr <- fish((pf_ch * caw_gr), t, 1) / sum(fish((pf_ch * caw_gr), t, 1))

fex_ch_p <- list(fex_ch_kn, fex_ch_br, fex_ch_cd, fex_ch_gr)
gex_ch_p <- list(gex_ch_kn, gex_ch_br, gex_ch_cd, gex_ch_gr)

# kuiper
# fines exit pmfs (slow fraction)
fex_kp_kn <- fish(1-(pf_kp * caw_kn), t, 1) / sum(fish(1-(pf_kp * caw_kn), t, 1))
fex_kp_br <- fish(1-(pf_kp * caw_br), t, 1) / sum(fish(1-(pf_kp * caw_br), t, 1))
fex_kp_cd <- fish(1-(pf_kp * caw_cd), t, 1) / sum(fish(1-(pf_kp * caw_cd), t, 1))
fex_kp_gr <- fish(1-(pf_kp * caw_gr), t, 1) / sum(fish(1-(pf_kp * caw_gr), t, 1))

# gravels (slow fraction)
gex_kp_kn <- fish((pf_kp * caw_kn), t, 1) / sum(fish((pf_kp * caw_kn), t, 1))
gex_kp_br <- fish((pf_kp * caw_br), t, 1) / sum(fish((pf_kp * caw_br), t, 1))
gex_kp_cd <- fish((pf_kp * caw_cd), t, 1) / sum(fish((pf_kp * caw_cd), t, 1))
gex_kp_gr <- fish((pf_kp * caw_gr), t, 1) / sum(fish((pf_kp * caw_gr), t, 1))

fex_kp_p <- list(fex_kp_kn, fex_kp_br, fex_kp_cd, fex_kp_gr)
gex_kp_p <- list(gex_kp_kn, gex_kp_br, gex_kp_cd, gex_kp_gr)

# kolmogorov-smirnov
# fines exit pmfs (slow fraction)
fex_ks_kn <- fish(1-(pf_ks * caw_kn), t, 1) / sum(fish(1-(pf_ks * caw_kn), t, 1))
fex_ks_br <- fish(1-(pf_ks * caw_br), t, 1) / sum(fish(1-(pf_ks * caw_br), t, 1))
fex_ks_cd <- fish(1-(pf_ks * caw_cd), t, 1) / sum(fish(1-(pf_ks * caw_cd), t, 1))
fex_ks_gr <- fish(1-(pf_ks * caw_gr), t, 1) / sum(fish(1-(pf_ks * caw_gr), t, 1))

# gravels (slow fraction)
gex_ks_kn <- fish((pf_ks * caw_kn), t, 1) / sum(fish((pf_ks * caw_kn), t, 1))
gex_ks_br <- fish((pf_ks * caw_br), t, 1) / sum(fish((pf_ks * caw_br), t, 1))
gex_ks_cd <- fish((pf_ks * caw_cd), t, 1) / sum(fish((pf_ks * caw_cd), t, 1))
gex_ks_gr <- fish((pf_ks * caw_gr), t, 1) / sum(fish((pf_ks * caw_gr), t, 1))

fex_ks_p <- list(fex_ks_kn, fex_ks_br, fex_ks_cd, fex_ks_gr)
gex_ks_p <- list(gex_ks_kn, gex_ks_br, gex_ks_cd, gex_ks_gr)






# bootstrap gravel ages 2
gravel_synth <- function(window = .25) {
  gad <- 0
  gch <- 0
  gkp <- 0
  gks <- 0

  for (i in 1:length(dfi)) {
    gad[i] <- sample(dfi, 1, prob = dfp)
    gch[i] <- sample(dfi, 1, prob = dfp)
    gkp[i] <- sample(dfi, 1, prob = dfp)
    gks[i] <- sample(dfi, 1, prob = dfp)
    rd <- runif(1)
    if (rd <= 0.5 - window) {

      gad[i] <- gad[i] - sample(tad, 1, prob = pexb_ad)
      gch[i] <- gch[i] - sample(tch, 1, prob = pexb_ch)
      gkp[i] <- gkp[i] - sample(tkp, 1, prob = pexb_kp)
      gks[i] <- gks[i] - sample(tks, 1, prob = pexb_ks)
    }
    if (rd > 0.5 + window) {
      pr <- sample(list(pex_kn, pex_br, pex_cd, pex_gr), 1, replace = T, prob = c(gct_kn, gct_br, gct_cd, gct_gr) / gct_tl)[[1]]
      pad <- sample(gex_ad_p, 1, replace = T, prob = gct_p)[[1]]
      pch <- sample(gex_ch_p, 1, replace = T, prob = gct_p)[[1]]
      pkp <- sample(gex_kp_p, 1, replace = T, prob = gct_p)[[1]]
      pks <- sample(gex_ks_p, 1, replace = T, prob = gct_p)[[1]]

      gad[i] <- gad[i] + sample(index, 1, replace = T, prob = pad)
      gch[i] <- gch[i] + sample(index, 1, replace = T, prob = pch)
      gkp[i] <- gkp[i] + sample(index, 1, replace = T, prob = pkp)
      gks[i] <- gks[i] + sample(index, 1, replace = T, prob = pks)
    }
  }
  list(gad, gch, gkp, gks)
}



gravel_fit <- function(n = 10, window = .25) {
  fits <- gravel_synth(window)
  gad <- fits[[1]]
  gch <- fits[[2]]
  gkp <- fits[[3]]
  gks <- fits[[4]]
  for (i in 2:n) {
    fits <- gravel_synth()
    gad <- c(gad, fits[[1]])
    gch <- c(gch, fits[[2]])
    gkp <- c(gkp, fits[[3]])
    gks <- c(gks, fits[[4]])
  }
  list(emp_cdf(gad), emp_cdf(gch), emp_cdf(gkp), emp_cdf(gks))
}

gf <- gravel_fit()
lines(gf[[1]])
lines(gf[[2]])
lines(gf[[3]])
lines(gf[[4]])



# bootstrap fines ages 2
fines_synth <- function(gf, # element from list output by gravel_synth (match to test)
                        prs, # probability list for slow fraction (match to test)
                        prf, # prob (fast fraction) (match to test)
                        prfi, # index for fast fraction probs
                        gcp = gct_p, # reach length selection probability
                        window = .25) {
  df <- charcoal$mn[charcoal$facies == 'DF'] # debris flow age obs
  dfc <- emp_cdf(df) # cdf of ages
  dfi <- dfc[ , 1] # index of cdf ages
  dfp <- to_pmf(dfc[ , 2]) # pmf of index ages
  s <- 0 # synth deposits
  for (i in 1:length(dfi)) {
    # pull from debris-flow or gravel distribution
    s[i] <- sample(sample(list(df, gf), 1)[[1]], 1)
    rd <- runif(1)
    if (rd <= 0.5 - window) { # fast fraction
      s[i] <- s[i] - sample(prfi, 1, prob = prf)
    }
    if (rd > 0.5 + window) { # slow fraction
      p <- sample(prs, 1, replace = T, prob = gcp)[[1]]
      s[i] <- s[i] + sample(index, 1, replace = T, prob = p)
    }
  }
  s
}

fines_fit <- function(gf, # element from list output by gravel_synth (match to test)
                      prs, # probability list for slow fraction (match to test)
                      prf, # prob (fast fraction) (match to test)
                      prfi, # index for fast fraction probs
                      gcp = gct_p, # reach length selection probability
                      n = 10, # number of times run simulation
                      window = .25) {
  s <- 0
  for (i in 1:n) {
    if (i == 1) {
      s <- fines_synth(gf, prs, prf, prfi, gcp, window)
    } else {
      s <- c(s, fines_synth(gf, prs, prf, prfi, gcp, window))
    }
  }
  emp_cdf(s)
}

gravel_synth <- function(window = .25) {
  df <- charcoal$mn[charcoal$facies == 'DF'] # debris flow age obs
  dfc <- emp_cdf(df) # cdf of ages
  dfi <- dfc[ , 1] # index of cdf ages
  dfp <- to_pmf(dfc[ , 2]) # pmf of index ages
  s <- 0 # synth deposits

  for (i in 1:length(dfi)) {
    s[i] <- sample(dfi, 1, prob = dfp)
    rd <- runif(1)
    if (rd <= 0.5 - window) {
      s[i] <- s[i] - sample(prfi, 1, prob = prf)
    }
    if (rd > 0.5 + window) {
      p <- sample(prs, 1, replace = T, prob = gcp)[[1]]
      s[i] <- s[i] + sample(index, 1, replace = T, prob = p)
    }
  }
  list(gad, gch, gkp, gks)
}

gravel_fit <- function(n = 10, window = .25) {
  fits <- gravel_synth(window)
  gad <- fits[[1]]
  gch <- fits[[2]]
  gkp <- fits[[3]]
  gks <- fits[[4]]
  for (i in 2:n) {
    fits <- gravel_synth()
    gad <- c(gad, fits[[1]])
    gch <- c(gch, fits[[2]])
    gkp <- c(gkp, fits[[3]])
    gks <- c(gks, fits[[4]])
  }
  list(emp_cdf(gad), emp_cdf(gch), emp_cdf(gkp), emp_cdf(gks))
}

# bootstrap gravels for plotting
gr_ad <- fines_fit(charcoal$mn[charcoal$facies == 'DF'], gex_ad_p, fexb_ad, tad, n = 1, window = .25)
gr_ch <- fines_fit(charcoal$mn[charcoal$facies == 'DF'], gex_ch_p, fexb_ch, tch, n = 1, window = .25)
gr_kp <- fines_fit(charcoal$mn[charcoal$facies == 'DF'], gex_kp_p, fexb_kp, tkp, n = 1, window = .25)
gr_ks <- fines_fit(charcoal$mn[charcoal$facies == 'DF'], gex_ks_p, fexb_ks, tks, n = 1, window = .25)

gr_fit <- gravel_fit(n = 10, win = fep^3)
gr_ad <- gr_fit[[1]]
gr_ch <- gr_fit[[2]]
gr_kp <- gr_fit[[3]]
gr_ks <- gr_fit[[4]]


png('gravel_fits.png', height = 17, width = 21, units = 'cm', res = 300)
magplot(emp_cdf(charcoal$mn[charcoal$facies == 'FF']), pch = 20, col = get_palette('coral'),
     xlab = 'Charcoal Age [y]', ylab = 'CDF', xlim = c(0, 17000))
points(dfc, pch = 20, col = get_palette('hardwood'))
points(emp_cdf(charcoal$mn[charcoal$facies == 'FG']), pch = 20, col = get_palette('charcoal'))
lines(gr_ad, lwd = 2.5, col = get_palette('ocean'))
lines(gr_ch, lwd = 2.5, col = get_palette('gold'))
lines(gr_kp, lwd = 2.5, col = get_palette('violet'))
lines(gr_ks, lwd = 2.5, col = get_palette('crimson'))
# lines(gf[[1]], lwd = 2.5, col = get_palette('ocean'))
# lines(gf[[2]], lwd = 2.5, col = get_palette('gold'))
# lines(gf[[3]], lwd = 2.5, col = get_palette('violet'))
# lines(gf[[4]], lwd = 2.5, col = get_palette('crimson'))
legend('bottomright', legend = c('Anderson-Darling', 'Chi-squared', 'Kuiper', 'Kolmogorov-Smirnov', 'Debris Flows', 'Gravels', 'Fines'),
       pch = c(rep(NA, 4), rep(20, 3)), lty = c(rep(1, 4), rep(NA, 3)),
       col = get_palette(c('ocean', 'gold', 'violet', 'crimson', 'hardwood', 'charcoal', 'coral'), .77))
dev.off()




# bootstrap fines for plotting
win = .12
fn_ad <- fines_fit(gf[[1]], fex_ad_p, fexb_ad, tad, n = 1, window = .12)
fn_ch <- fines_fit(gf[[2]], fex_ch_p, fexb_ch, tch, n = 1, window = .12)
fn_kp <- fines_fit(gf[[3]], fex_kp_p, fexb_kp, tkp, n = 1, window = .09)
fn_ks <- fines_fit(gf[[4]], fex_ks_p, fexb_ks, tks, n = 1, window = .09)

fn_fit <- fine_fit(mod = ffep^2)
fn_ad <- fn_fit[[1]]
fn_ch <- fn_fit[[2]]
fn_kp <- fn_fit[[3]]
fn_ks <- fn_fit[[4]]


png('fines_fits.png', height = 17, width = 21, units = 'cm', res = 300)
magplot(emp_cdf(charcoal$mn[charcoal$facies == 'FF']), pch = 20, col = get_palette('coral'),
     xlab = 'charcoal age', ylab = 'CDF of samples', xlim = c(0, 17000))
points(dfc, pch = 20, col = get_palette('hardwood'))
points(emp_cdf(charcoal$mn[charcoal$facies == 'FG']), pch = 20, col = get_palette('charcoal'))
lines(fn_ad, lwd = 2.5, col = get_palette('ocean'))
lines(fn_ch, lwd = 2.5, col = get_palette('gold'))
lines(fn_kp, lwd = 2.5, col = get_palette('violet'))
lines(fn_ks, lwd = 2.5, col = get_palette('crimson'))
legend('bottomright', legend = c('anderson-darling', 'chi-squared', 'kuiper', 'kolmogorov-smirnov', 'debris flows', 'gravels', 'fines'),
       pch = c(rep(NA, 4), rep(20, 3)), lty = c(rep(1, 4), rep(NA, 3)),
       col = get_palette(c('ocean', 'gold', 'violet', 'crimson', 'hardwood', 'charcoal', 'coral'), .77))
dev.off()


gvl <- data.table::fread('/home/erik/output/gravels_test.csv')
gvl <- unlist(gvl)
magplot(emp_cdf(charcoal$mn[charcoal$facies == 'FF']), pch = 20, col = get_palette('coral'),
        xlab = 'charcoal age', ylab = 'CDF of samples', xlim = c(0, 17000))
points(dfc, pch = 20, col = get_palette('hardwood'))
points(emp_cdf(charcoal$mn[charcoal$facies == 'FG']), pch = 20, col = get_palette('charcoal'))
lines(emp_cdf(gvl))



hplot(emp_cdf(pebs$b), log = 'x', xlim = c(4.5, 7.5))

b4m$coefficients[1] + (b4m$coefficients[2] * log(7.4))




tc <- dfi[dfc[, 2] > (1 - te / (te + tr))]
df_tc <- dfp[dfc[, 2] > (1 - te / (te + tr))]
cpex <- pex_pmf[index >= min(tc)]
cpexi <- index[index >= min(tc)]
cfl <- rep(0, length(tc))
cfl[1] <- sum(cpex[cpexi < tc[1]])
for (i in 2:length(tc)) {
  cfl[i] <- sum(cpex[cpexi > tc[i-1] & cpexi <= tc[i]])
}
cfl <- c(rep(0,length(dfi)-length(tc)), cfl)
cfl <- convo_add(cfl, df_tc, tc)
tfl <- c(dfp[dfc[, 2] <= (1 - te / (te + tr))], cfl)
plot(dfi, cumsum(tfl))
plot(tc, cumsum(cfl))

cfl <- rep(0, length(dfi))
cfl[1] <- sum(pex_pmf[index < dfi[1]])
for (i in 2:length(dfi)) {
  cfl[i] <- sum(pex_pmf[index > dfi[i-1] & index <= tc[i]])
}
# cfl <- convo_add(cfl * (tr/(te+tr)), dfp, dfi)
cfl <- convo_add(dfp, cfl, dfi)
cfl[is.na(cfl)] <- 0
cfl <- cfl / sum(cfl)
tfl <- cfl * 0.5 * tr/(te+tr) + dfp * te/(te+tr)
tfl <- tfl / sum(tfl)
plot(dfi, cumsum(cfl))
plot(dfi, cumsum(tfl))

lines(emp_cdf(charcoal$mn[charcoal$facies == 'DF']))
lines(emp_cdf(charcoal$mn[charcoal$facies == 'FG']))

cflb <- rep(0, length(dfi))
cflb[1] <- sum(pex_b[t * tp < dfi[1]])
for (i in 2:length(dfi)) {
  cflb[i] <- sum(pex_b[(t * tp) > dfi[i-1] & (t * tp) <= dfi[i]])
}

cflb <- convo(cfl, cflb, dfi)
cflb <- cflb / sum(cflb)
tfl <- convo_add(cfl,  cflb, dfi)
tfl <- tfl / sum(tfl)
plot(dfi, cumsum(tfl))

plot(index, pex_pmf)
plot(dfi, dfp, log = 'y')
plot(dfc)
lines(index, pex_pmf)

lines(index, cumsum(df_pmf))
plot(dfi, cumsum(tfl))
lines(dfc)
lines(emp_cdf(charcoal$mn[charcoal$facies == 'FG']))
plot(dfi, cumsum(cflb))
plot(index, cumsum(df_pmf))
plot(index, cumsum(pex_b))
lines(index, cumsum(df_pmf))
lines(dfi, cumsum(dfp), col = 'red', lwd = 3)
lines(emp_cdf(charcoal$mn[charcoal$facies == 'DF']))



png('predicted_fluvial_age.png', height = 20, width = 24, units = 'cm', res = 300)
plot(dfi, tfl %>% to_cdf, type = 'l', col = get_palette('crimson'), lwd = 3,
     xlab = 'charcoal age [years]',
     ylab = 'CDF', log = 'x', ylim = c(0,1))
lines(emp_cdf(charcoal$mn[charcoal$facies == 'FF']), col = get_palette('sky'), lwd = 3)
lines(emp_cdf(charcoal$mn[charcoal$facies == 'FG']), col = get_palette('ocean'), lwd = 3)
lines(emp_cdf(charcoal$mn[charcoal$facies == 'DF']), col = get_palette('hardwood'), lwd = 3)
lines(emp_cdf(df_stereo), col = get_palette('crimson'), lwd = 3)
legend('bottomright', legend = c('gravels', 'fines', 'predicted'),
       fill = get_palette(c('ocean', 'sky', 'crimson'), .8))
dev.off()


# average traversal distance, including tribs
knsa <- sf::st_read("/media/erik/catacomb/gis/creeks/knowles_sa.shp")
brsa <- sf::st_read("/media/erik/catacomb/gis/creeks/bear_sa.shp")
plot(brsa)

# add main stem to travel distance of tribs
# from headwater
# k1 node id: 384365-384391, junction id: 384327
knsa$ToMouth_km[knsa$NODE_ID %in% 384365:384391] <- knsa$ToMouth_km[knsa$NODE_ID %in% 384365:384391] + knsa$ToMouth_km[knsa$NODE_ID == 384327]

# k2 node id: 3843392-384395, junction id: 384228
knsa$ToMouth_km[knsa$NODE_ID %in% 384392:384395] <- knsa$ToMouth_km[knsa$NODE_ID %in% 384392:384395] + knsa$ToMouth_km[knsa$NODE_ID == 384228]

# k3 node id: 3843396-384405, junction id: 384214
knsa$ToMouth_km[knsa$NODE_ID %in% 384396:384405] <- knsa$ToMouth_km[knsa$NODE_ID %in% 384396:384405] + knsa$ToMouth_km[knsa$NODE_ID == 384214]

# k4 node id: 384406-384481, junction id: 384190
knsa$ToMouth_km[knsa$NODE_ID %in% 384406:384481] <- knsa$ToMouth_km[knsa$NODE_ID %in% 384406:384481] + knsa$ToMouth_km[knsa$NODE_ID == 384190]

# k4a node id: 384482-384509, junction id: 384454
knsa$ToMouth_km[knsa$NODE_ID %in% 384482:384509] <- knsa$ToMouth_km[knsa$NODE_ID %in% 384482:384509] + knsa$ToMouth_km[knsa$NODE_ID == 384454]

# k4b node id: 384510-384527, junction id: 384430
knsa$ToMouth_km[knsa$NODE_ID %in% 384510:384527] <- knsa$ToMouth_km[knsa$NODE_ID %in% 384510:384527] + knsa$ToMouth_km[knsa$NODE_ID == 384430]

# k5 node id: 384528:384536, junction id: 384162
knsa$ToMouth_km[knsa$NODE_ID %in% 384528:384536] <- knsa$ToMouth_km[knsa$NODE_ID %in% 384528:384536] + knsa$ToMouth_km[knsa$NODE_ID == 384162]

# k6 node id: 384537:384657, junction id: 384145
knsa$ToMouth_km[knsa$NODE_ID %in% 384537:384657] <- knsa$ToMouth_km[knsa$NODE_ID %in% 384537:384657] + knsa$ToMouth_km[knsa$NODE_ID == 384145]

# k6a node id: 384658:384668, junction id: 384647
knsa$ToMouth_km[knsa$NODE_ID %in% 384658:384668] <- knsa$ToMouth_km[knsa$NODE_ID %in% 384658:384668] + knsa$ToMouth_km[knsa$NODE_ID == 384647]

# k6b node id: 384669:384684, junction id: 384635
knsa$ToMouth_km[knsa$NODE_ID %in% 384669:384684] <- knsa$ToMouth_km[knsa$NODE_ID %in% 384669:384684] + knsa$ToMouth_km[knsa$NODE_ID == 384635]

# k6c node id: 384685:384721, junction id: 384615
knsa$ToMouth_km[knsa$NODE_ID %in% 384685:384721] <- knsa$ToMouth_km[knsa$NODE_ID %in% 384685:384721] + knsa$ToMouth_km[knsa$NODE_ID == 384615]

# k6d node id: 384722:384738, junction id: 384614
knsa$ToMouth_km[knsa$NODE_ID %in% 384722:384738] <- knsa$ToMouth_km[knsa$NODE_ID %in% 384722:384738] + knsa$ToMouth_km[knsa$NODE_ID == 384614]

# k6e node id: 384739:384758, junction id: 384577
knsa$ToMouth_km[knsa$NODE_ID %in% 384739:384758] <- knsa$ToMouth_km[knsa$NODE_ID %in% 384739:384758] + knsa$ToMouth_km[knsa$NODE_ID == 384577]

# k7 node id: 384759:384782, junction id: 384111
knsa$ToMouth_km[knsa$NODE_ID %in% 384759:384782] <- knsa$ToMouth_km[knsa$NODE_ID %in% 384759:384782] + knsa$ToMouth_km[knsa$NODE_ID == 384111]

# k8 node id: 384783:384806, junction id: 384111
knsa$ToMouth_km[knsa$NODE_ID %in% 384783:384806] <- knsa$ToMouth_km[knsa$NODE_ID %in% 384783:384806] + knsa$ToMouth_km[knsa$NODE_ID == 384111]

# k8 node id: 384807:384941, junction id: 384041
knsa$ToMouth_km[knsa$NODE_ID %in% 384807:384941] <- knsa$ToMouth_km[knsa$NODE_ID %in% 384807:384941] + knsa$ToMouth_km[knsa$NODE_ID == 384041]

# k8a node id: 384942:384947, junction id: 384871
knsa$ToMouth_km[knsa$NODE_ID %in% 384942:384947] <- knsa$ToMouth_km[knsa$NODE_ID %in% 384942:384947] + knsa$ToMouth_km[knsa$NODE_ID == 384871]

# k9 node id: 384948:384964, junction id: 384024
knsa$ToMouth_km[knsa$NODE_ID %in% 384948:384964] <- knsa$ToMouth_km[knsa$NODE_ID %in% 384948:384964] + knsa$ToMouth_km[knsa$NODE_ID == 384024]

# k10 node id: 384965:385004, junction id: 384009
knsa$ToMouth_km[knsa$NODE_ID %in% 384965:385004] <- knsa$ToMouth_km[knsa$NODE_ID %in% 384965:385004] + knsa$ToMouth_km[knsa$NODE_ID == 384009]

# k11 node id: 385005:385070, junction id: 384004
knsa$ToMouth_km[knsa$NODE_ID %in% 385005:385070] <- knsa$ToMouth_km[knsa$NODE_ID %in% 385005:385070] + knsa$ToMouth_km[knsa$NODE_ID == 384004]

# k11a node id: 385071:385082, junction id: 385045
knsa$ToMouth_km[knsa$NODE_ID %in% 385071:385082] <- knsa$ToMouth_km[knsa$NODE_ID %in% 385071:385082] + knsa$ToMouth_km[knsa$NODE_ID == 385045]

# adjust starting distance to beginning of study area
knsa <- knsa[knsa$NODE_ID >= min(creeks$NODE_ID[creeks$creek_name == 'knowles']) &
               knsa$NODE_ID < 385082, ]
knsa$dist <- knsa$ToMouth_km - min(knsa$ToMouth_km)

# travel distance for bear creek
# b1 node id: 386453:386459, junction id: 386440
brsa$ToMouth_km[brsa$NODE_ID %in% 386453:386459] <- brsa$ToMouth_km[brsa$NODE_ID %in% 386453:386459] + brsa$ToMouth_km[brsa$NODE_ID == 386440]

# b2 node id: 386460:386489, junction id: 386409
brsa$ToMouth_km[brsa$NODE_ID %in% 386460:386489] <- brsa$ToMouth_km[brsa$NODE_ID %in% 386460:386489] + brsa$ToMouth_km[brsa$NODE_ID == 386409]

# b3 node id: 386490:386503, junction id: 386389
brsa$ToMouth_km[brsa$NODE_ID %in% 386490:386503] <- brsa$ToMouth_km[brsa$NODE_ID %in% 386490:386503] + brsa$ToMouth_km[brsa$NODE_ID == 386389]

# b4 node id: 386504:386553, junction id: 386378
brsa$ToMouth_km[brsa$NODE_ID %in% 386504:386553] <- brsa$ToMouth_km[brsa$NODE_ID %in% 386504:386553] + brsa$ToMouth_km[brsa$NODE_ID == 386378]

# b5 node id: 386554:386610, junction id: 386359
brsa$ToMouth_km[brsa$NODE_ID %in% 386554:386610] <- brsa$ToMouth_km[brsa$NODE_ID %in% 386554:386610] + brsa$ToMouth_km[brsa$NODE_ID == 386359]

# b5a node id: 386611:386645, junction id: 386574
brsa$ToMouth_km[brsa$NODE_ID %in% 386611:386645] <- brsa$ToMouth_km[brsa$NODE_ID %in% 386611:386645] + brsa$ToMouth_km[brsa$NODE_ID == 386574]

# b6 node id: 386646:386659, junction id: 386348
brsa$ToMouth_km[brsa$NODE_ID %in% 386646:386659] <- brsa$ToMouth_km[brsa$NODE_ID %in% 386646:386659] + brsa$ToMouth_km[brsa$NODE_ID == 386348]

# b7 node id: 386660:386693, junction id: 386331
brsa$ToMouth_km[brsa$NODE_ID %in% 386660:386693] <- brsa$ToMouth_km[brsa$NODE_ID %in% 386660:386693] + brsa$ToMouth_km[brsa$NODE_ID == 386331]

# b8 node id: 386694:386713, junction id: 386314
brsa$ToMouth_km[brsa$NODE_ID %in% 386694:386713] <- brsa$ToMouth_km[brsa$NODE_ID %in% 386694:386713] + brsa$ToMouth_km[brsa$NODE_ID == 386314]

# b9 node id: 386714:386727, junction id: 386299
brsa$ToMouth_km[brsa$NODE_ID %in% 386714:386727] <- brsa$ToMouth_km[brsa$NODE_ID %in% 386714:386727] + brsa$ToMouth_km[brsa$NODE_ID == 386299]

# b10 node id: 386728:386744, junction id: 386278
brsa$ToMouth_km[brsa$NODE_ID %in% 386728:386744] <- brsa$ToMouth_km[brsa$NODE_ID %in% 386728:386744] + brsa$ToMouth_km[brsa$NODE_ID == 386278]

# b11 node id: 386745:386762, junction id: 386269
brsa$ToMouth_km[brsa$NODE_ID %in% 386745:386762] <- brsa$ToMouth_km[brsa$NODE_ID %in% 386745:386762] + brsa$ToMouth_km[brsa$NODE_ID == 386269]

# b12 node id: 386763:386777, junction id: 386258
brsa$ToMouth_km[brsa$NODE_ID %in% 386763:386777] <- brsa$ToMouth_km[brsa$NODE_ID %in% 386763:386777] + brsa$ToMouth_km[brsa$NODE_ID == 386258]

# b13 node id: 386778:386796, junction id: 386243
brsa$ToMouth_km[brsa$NODE_ID %in% 386778:386796] <- brsa$ToMouth_km[brsa$NODE_ID %in% 386778:386796] + brsa$ToMouth_km[brsa$NODE_ID == 386243]

# b14 node id: 386797:386814, junction id: 386239
brsa$ToMouth_km[brsa$NODE_ID %in% 386797:386814] <- brsa$ToMouth_km[brsa$NODE_ID %in% 386797:386814] + brsa$ToMouth_km[brsa$NODE_ID == 386239]

plot(emp_cdf(knsa$dist))
setwd('/media/erik/catacomb/research/')
png('travel_dist.png', height = 17, width = 21, units = 'cm', res = 300)
plot(density(knsa$dist),
     xlab = 'Travel Distance [km]',
     main = '', ylim = c(0, .7), col = get_palette('crimson', .7), lwd = 2.5)
lines(density(brsa$ToMouth_km), col = get_palette('ocean', 0.7), lwd = 2.5)
lines(density(c(brsa$ToMouth_km, knsa$dist)), col = get_palette('charcoal', 0.7), lwd = 2.5)
abline(v = mean(c(brsa$ToMouth_km, knsa$dist)), col = get_palette('charcoal', .7), lty = 2, lwd = 1.5)
abline(v = mean(c(brsa$ToMouth_km)), lty = 2, col = get_palette('ocean', .7), lwd = 1.5)
abline(v = mean(c(knsa$dist)), lty = 2, col = get_palette('crimson', .7), lwd = 1.5)
legend('topright', legend = c('Bear Creek', 'Knowles Creek', 'Both Creeks', 'Mean'),
       col = get_palette(c('ocean', 'crimson', 'charcoal', 'charcoal'), 0.7), lty = c(1, 1, 1, 2))
dev.off()

mn_d <- mean(c(brsa$ToMouth_km, knsa$dist))
mn_dkn <- mean(knsa$dist)


pcdf <- emp_cdf(pebs$b)
d90s <- pcdf[pcdf[,2] > .9,1] %>% min
d90s_mass <- exp(predict(gmod, newdata = data.frame(b = d90s)))
st_d90 <- -(log(d2_mass) - log(d90s_mass)) / mn_dkn
st_d90s <- -(log(d2_mass) - log(d50s_mass)) / mn_dkn
st_dst <- -(log(d2_mass) - log(.1)) / mn_dkn

plot(emp_cdf(surf))
plot(emp_cdf(exp(predict(gmod, newdata = data.frame(b = surf)))))
plot(emp_cdf(exp(predict(gmod, newdata = data.frame(b = surf)) * exp(-stm * mn_dkn))))
mass_rem <- -(log(d2_mass) - predict(gmod, newdata = data.frame(b = surf))) - stm * mn_dkn
mass_rem[mass_rem < 0] <- 0
sum(mass_rem)
sum(exp(predict(gmod, newdata = data.frame(b = surf))))
plot(emp_cdf(exp(predict(gmod, newdata = data.frame(b = surf)))))
lines(emp_cdf(exp(mass_rem)))
plot(emp_cdf(predict(gmod, newdata = data.frame(b = surf)) - mass_rem))

st_st <- pcdf[pcdf[,1] <= 4.69, ]
st_st[st_st[,1] > 4.68, 2]

df_trans_ks <- data.table::fread('/home/erik/output/transits_cdf_gr_kp.csv')
gr_trans_ks <- data.table::fread('/home/erik/output/transits_cdf_gr_ks.csv')
gr_trans_ks <- data.table::fread('/home/erik/output/gravels_transits_ks.csv')
gr_trans_ks <- data.table::fread('/home/erik/output/gravels_transits_kp.csv')
gr_trans_kp <- data.table::fread('/home/erik/output/transits_cdf_gr_kp.csv')
gr_trans_ch <- data.table::fread('/home/erik/output/transits_cdf_gr_ch.csv')
gr_trans_ks <- data.table::fread('/home/erik/output/transits_cdf_fn_ks.csv')
gr_trans_kp <- data.table::fread('/home/erik/output/transits_cdf_fn_kp.csv')
df_trans_ks <- data.table::fread('/home/erik/output/debris_flow_deposits_cdf.csv')


png('gravels_fit.png', height = 17, width = 21, units = 'cm', res = 300)
magicaxis::magplot(emp_cdf(charcoal$mn[charcoal$facies == 'FF'] + 50), pch = 20, col = get_palette('coral'),
                   xlab = 'Charcoal Age [y]', ylab = 'CDF of samples', xlim = c(0, 17000))
points(emp_cdf(debris_flows), pch = 20, col = get_palette('hardwood'))
points(emp_cdf(charcoal$mn[charcoal$facies == 'FG'] + 50), pch = 20, col = get_palette('charcoal'))

lines(0:20000, cumsum(unlist(gr_trans_ks)), lwd = 2.5, col = get_palette('crimson', .7))
lines(0:20000, cumsum(unlist(gr_trans_kp)), lwd = 2.5, col = get_palette('violet', .7))
# lines(0:20000, cumsum(unlist(gr_trans_ch)), lwd = 2.5, col = get_palette('gold', .7))
legend('bottomright', legend = c('Debris-Flows', 'Gravels', 'Fines', 'Kolmogorov-Smirnov', 'Kuiper'),
       lty = c(NA, NA, NA, 1, 1), lwd = c(NA, NA, NA, 2, 2),
       pch = c(20, 20, 20, NA, NA), col = get_palette(c('hardwood', 'charcoal', 'coral', 'crimson', 'violet'), .9))
dev.off()


head(df_trans_ks)
df_trans_ks[df_trans_ks < 0] <- 0
write.csv(df_trans_ks, 'gravels_transits_kp.csv', row.names = F)



png('pebble_counts.png', height = 17, width = 21, units = 'cm', res = 300)
magicaxis::magplot(emp_cdf(pebs$b), log = 'x', pch = 20,
                   col = get_palette('ocean'), xlim = c(min(pebs$b), max(surf)),
                   xlab = 'b-axis width [mm]',
                   ylab = 'CDF of samples')
points(emp_cdf(surf), pch = 20, col = get_palette('crimson'))
abline(v = median(pebs$b), lty = 2)
text(9.5, 0.45, labels = paste0('D50s = ', round(median(pebs$b), 2), ' g'))
abline(v = median(surf), lty = 2)
text(65, 0.45, labels = paste0('D50 = ', round(median(surf), 2), ' g'))
legend('bottomright', legend = c('Surface', 'Subsurface', 'Median'),
       pch = c(20, 20, NA), lty = c(NA, NA, 2), lwd = c(NA, NA, 1),
       col = get_palette(c('crimson', 'ocean', 'charcoal'), .9))
dev.off()


