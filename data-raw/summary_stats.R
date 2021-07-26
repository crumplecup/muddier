
transits <- data.table::fread('/home/erik/output/transit_fits.csv')
plot(transits$rate, transits$mean)


slp <- log(creeks$slope[creeks$slope > 0 &
                          creeks$creek_name != 'cedar' &
                          creeks$creek_name != 'hoffman'])
ca <- log(creeks$contr_area[creeks$slope > 0 &
                              creeks$creek_name != 'cedar' &
                              creeks$creek_name != 'hoffman'])

# signs of data quality issues in cedar
# ced <- creeks[creeks$creek_name == 'cedar', ]
# ced <- ced[1:111, ]
# slp <- c(slp, log(ced$slope))
# ca <- c(ca, log(ced$contr_area))

y <- lm(slp ~ ca)
summary(y)
plot(slp, ca)

# plot log slope vs contr area for each creek
# scale invariance evident in knowles
kn <- creeks[creeks$creek_name == 'knowles' & creeks$slope > 0, ]
kn_fit <- lm(log(contr_area) ~ log(slope), data = kn)
kn_pred <- predict(kn_fit, newdata = kn)
kn_epred <- exp(kn_pred)
kn_mod <- lm(log(slope) ~ log(contr_area), data = kn)
kn_s <- summary(kn_mod)
setwd('/media/erik/catacomb/research/')
png('kn_slope_to_area.png', height = 16, width = 20, units = 'cm', res = 300)
plot(kn$slope, kn$contr_area, log = 'xy',
     xlab = 'slope', ylab = 'contributing area (km2)',
     pch = 20, col = get_palette('ocean'))
lines(kn$slope, kn_epred, lwd = 3, col = get_palette('charcoal', .7))
text(0.12, 1.2, paste0('R2 = ' , round(kn_s$adj.r.squared, 3)))
text(0.12, 1.0, paste0('β = ' , round(kn_s$coefficients[2, 1], 3)))
dev.off()

# bear
br <- creeks[creeks$creek_name == 'bear' & creeks$slope > 0, ]
br_fit <- lm(log(contr_area) ~ log(slope), data = br)
br_pred <- predict(br_fit, newdata = br)
br_epred <- exp(br_pred)
br_mod <- lm(log(slope) ~ log(contr_area), data = br)
br_s <- summary(br_mod)
setwd('/media/erik/catacomb/research/')
png('br_slope_to_area.png', height = 16, width = 20, units = 'cm', res = 300)
plot(br$slope, br$contr_area, log = 'xy',
     xlab = 'slope', ylab = 'contributing area (km2)',
     pch = 20, col = get_palette('ocean'))
lines(br$slope, br_epred, lwd = 3, col = get_palette('charcoal', .7))
text(0.35, .25, paste0('R2 = ' , round(br_s$adj.r.squared, 3)))
text(0.35, .2, paste0('β = ' , round(br_s$coefficients[2, 1], 3)))
dev.off()

# cedar
cd <- creeks[creeks$creek_name == 'cedar' & creeks$slope > 0, ]
cd <- cd[1:111, ]
cd_fit <- lm(log(contr_area) ~ log(slope), data = cd)
cd_pred <- predict(cd_fit, newdata = cd)
cd_epred <- exp(cd_pred)
cd_mod <- lm(log(slope) ~ log(contr_area), data = cd)
cd_s <- summary(cd_mod)
setwd('/media/erik/catacomb/research/')
png('cd_slope_to_area.png', height = 16, width = 20, units = 'cm', res = 300)
plot(cd$slope, cd$contr_area, log = 'xy',
     xlab = 'slope', ylab = 'contributing area (km2)',
     pch = 20, col = get_palette('ocean'))
lines(cd$slope, cd_epred, lwd = 3, col = get_palette('charcoal', .7))
text(0.15, .7, paste0('R2 = ' , round(cd_s$adj.r.squared, 3)))
text(0.15, .6, paste0('β = ' , round(cd_s$coefficients[2, 1], 3)))
dev.off()

# hoffman
hf <- creeks[creeks$creek_name == 'hoffman', ]
hf <- hf[!is.na(hf$slope), ]
hf_fit <- lm(log(contr_area) ~ log(slope), data = hf)
hf_pred <- predict(hf_fit, newdata = hf)
hf_epred <- exp(hf_pred)
hf_mod <- lm(log(slope) ~ log(contr_area), data = hf)
hf_s <- summary(hf_mod)
setwd('/media/erik/catacomb/research/')
png('hf_slope_to_area.png', height = 16, width = 20, units = 'cm', res = 300)
plot(hf$slope, hf$contr_area, log = 'xy',
     xlab = 'slope', ylab = 'contributing area (km2)',
     pch = 20, col = get_palette('ocean'))
lines(hf$slope, hf_epred, lwd = 3, col = get_palette('charcoal', .7))
text(0.17, .2, paste0('R2 = ' , round(hf_s$adj.r.squared, 3)))
text(0.17, .16, paste0('β = ' , round(hf_s$coefficients[2, 1], 3)))
dev.off()

# creeks
cr <- creeks[!is.na(creeks$slope) & creeks$slope > 0 & creeks$creek_name != 'cedar', ]
cr_fit <- lm(log(contr_area) ~ log(slope), data = cr)
cr_pred <- predict(cr_fit, newdata = cr)
cr_epred <- exp(cr_pred)
cr_mod <- lm(log(slope) ~ log(contr_area), data = cr)
cr_s <- summary(cr_mod)
setwd('/media/erik/catacomb/research/')
png('cr_slope_to_area.png', height = 20, width = 24, units = 'cm', res = 300)
plot(cr$slope, cr$contr_area, log = 'xy',
     xlab = 'slope', ylab = 'contributing area (km2)',
     pch = 20, col = get_palette('slate', .01))
lines(kn$slope, kn_epred, lwd = 3, col = get_palette('ocean', .7))
points(kn$slope, kn$contr_area, pch = 20, col = get_palette('ocean'))
lines(br$slope, br_epred, lwd = 3, col = get_palette('hardwood', .7))
points(br$slope, br$contr_area, pch = 20, col = get_palette('hardwood'))
lines(hf$slope, hf_epred, lwd = 3, col = get_palette('crimson', .7))
points(hf$slope, hf$contr_area, pch = 20, col = get_palette('crimson'))
lines(cd$slope, cd_epred, lwd = 3, col = get_palette('forest', .7))
points(cd$slope, cd$contr_area, pch = 20, col = get_palette('forest', .5))
text(0.33, .25, paste0(
  'β = ' , round(br_s$coefficients[2, 1], 2), '  R2 = ' , round(br_s$adj.r.squared, 2)
  ), col = get_palette('hardwood', 1))
text(0.33, 1.0, paste0(
  'β = ' , round(cd_s$coefficients[2, 1], 2), '  R2 = ' , round(cd_s$adj.r.squared, 2)
), col = get_palette('forest', 1))
text(0.33, 0.055, paste0(
  'β = ' , round(hf_s$coefficients[2, 1], 2), '  R2 = ' , round(hf_s$adj.r.squared, 2)
), col = get_palette('crimson', 1))
text(0.33, 0.6, paste0(
  'β = ' , round(kn_s$coefficients[2, 1], 2), '  R2 = ' , round(kn_s$adj.r.squared, 2)
), col = get_palette('ocean', 1))
legend('bottomleft', legend = c('bear', 'cedar', 'hoffman', 'knowles'),
       fill = get_palette(c('hardwood', 'forest', 'crimson', 'ocean'), .9))
dev.off()

# plot log xsec area vs frequency for each creek
# scale invariance evident in knowles
cr <- creeks[creeks$xsec_area > 0, ]
# cdf of xsec area
fr <- emp_cdf(cr$xsec_area)
# add to creek object
cr$xsec_fr <- unlist(lapply(cr$xsec_area, function(x) max(fr[fr[,1] <= x, 2])))
# fit model
cr_fit <- lm(log(xsec_fr) ~ log(xsec_area), data = cr)
cr_s <- summary(cr_fit)
cr_mod <- lm(log(xsec_area) ~ log(xsec_fr), data = cr)
cr$pred <- exp(predict(cr_fit, newdata = cr))

kn <- cr[cr$creek_name == 'knowles', ]
kn <- kn[order(kn$xsec_area), ]
kn_fit <- lm(log(xsec_fr) ~ log(xsec_area), data = kn)
kn_s <- summary(kn_fit)
kn$pred <- exp(predict(kn_fit, newdata = kn))
br <- cr[cr$creek_name == 'bear', ]
br <- br[order(br$xsec_area), ]
br_fit <- lm(log(xsec_fr) ~ log(xsec_area), data = br)
br_s <- summary(br_fit)
br$pred <- exp(predict(br_fit, newdata = br))
hf <- cr[hf$creek_name == 'hoffman', ]
hf <- hf[order(hf$xsec_area), ]
hf_fit <- lm(log(xsec_fr) ~ log(xsec_area), data = hf)
hf_s <- summary(hf_fit)
hf$pred <- exp(predict(hf_fit, newdata = hf))
cd <- cr[cd$creek_name == 'cedar', ]
cd <- cd[order(cd$xsec_area), ]
cd_fit <- lm(log(xsec_fr) ~ log(xsec_area), data = cd)
cd_s <- summary(cd_fit)
cd$pred <- exp(predict(cd_fit, newdata = cd))

setwd('/media/erik/catacomb/research/')
png('cr_xsec_fr.png', height = 20, width = 24, units = 'cm', res = 300)
plot(cr$xsec_area, cr$xsec_fr, log = 'xy',
     xlab = 'cross-sectional area (m2)', ylab = 'CDF',
     pch = 20, col = get_palette('slate', .01), ylim = c(0.01, 1))
lines(kn$xsec_area, kn$pred, lwd = 3, col = get_palette('ocean', .3))
points(kn$xsec_area, kn$xsec_fr, pch = 20, col = get_palette('ocean', .1))
lines(br$xsec_area, br$pred, lwd = 3, col = get_palette('hardwood', .3))
points(br$xsec_area, br$xsec_fr, pch = 20, col = get_palette('hardwood', .1))
lines(hf$xsec_area, hf$pred, lwd = 3, col = get_palette('crimson', .3))
points(hf$xsec_area, hf$xsec_fr, pch = 20, col = get_palette('crimson', .035))
lines(cd$xsec_area, cd$pred, lwd = 3, col = get_palette('forest', .3))
points(cd$xsec_area, cd$xsec_fr, pch = 20, col = get_palette('forest', .03))
lines(cd$xsec_area, cd$pred, lwd = 3, col = get_palette('charcoal', .3))
text(0.2, 0.025, paste0(
  'β = ' , round(cr_s$coefficients[2, 1], 2), '  R2 = ' , round(cr_s$adj.r.squared, 2)
), col = get_palette('charcoal', 1))
text(0.2, 0.033, paste0(
  'β = ' , round(kn_s$coefficients[2, 1], 2), '  R2 = ' , round(kn_s$adj.r.squared, 2)
), col = get_palette('ocean', 1))
text(0.2, 0.043, paste0(
  'β = ' , round(hf_s$coefficients[2, 1], 2), '  R2 = ' , round(hf_s$adj.r.squared, 2)
), col = get_palette('crimson', 1))
text(0.2, 0.057, paste0(
  'β = ' , round(cd_s$coefficients[2, 1], 2), '  R2 = ' , round(cd_s$adj.r.squared, 2)
), col = get_palette('forest', 1))
text(0.2, 0.075, paste0(
  'β = ' , round(br_s$coefficients[2, 1], 2), '  R2 = ' , round(br_s$adj.r.squared, 2)
), col = get_palette('hardwood', 1))
legend('topleft', legend = c('bear', 'cedar', 'hoffman', 'knowles', 'all'),
       fill = get_palette(c('hardwood', 'forest', 'crimson', 'ocean', 'charcoal'), .9))
dev.off()

vfill <- sum(creeks$xsec_area * 10)  # m3
vfillh <- sum(creeks$xsec_area * 12.5)
ein <- max(creeks$contr_area[creeks$creek_name == 'bear']) +
  max(creeks$contr_area[creeks$creek_name == 'cedar']) +
  max(creeks$contr_area[creeks$creek_name == 'hoffman']) +
  max(creeks$contr_area[creeks$creek_name == 'knowles'])

sedvol <- data.table::fread('/home/erik/data/sedvol.csv')
sfill <- sum(sedvol$`Volume increment (m3)`[sedvol$Reach == 'Cedar']) +
  sum(sedvol$`Volume increment (m3)`[sedvol$Reach == 'Hoffman']) +
  sum(sedvol$`Volume increment (m3)`[sedvol$Reach %in% c('LB', 'UB')]) +
  sum(sedvol$`Volume increment (m3)`[sedvol$Reach %in% c('DFK', 'LK', 'UK')])
vfill - sfill

rec <- data.table::fread('/home/erik/output/df20k_1000.csv')
rec1 <- data.table::fread('/home/erik/output/df20k_1001.csv')
rec <- rbind(rec, rec1)
rm(rec1)
recf <- data.table::fread('/home/erik/output/ff20k_1000.csv')
recg <- data.table::fread('/home/erik/output/fg20k_1000.csv')


# Debris-flow transit time quantiles are
# [66.36173734474495, 301.13675474695447, 645.3447740909093, 1420.3701377105074, 6931.203620019353]
# Fluvial fines transit time quantiles are
# [67.80999012478165, 314.512768367902, 662.7312598279211, 1449.8351653715888, 6951.672643844576]
# Fluvial gravels transit time quantiles are
# [67.47723208784629, 319.6997800621376, 672.5117861952269, 1465.5612335698345, 6931.420597779025]

# debris-flow input/output rate
drt <- rec$input[rec$ks == min(rec$ks)[1]]
# fines input/output rate
frt <- recf$input[recf$ks == min(recf$ks)[1]][1]
# gravels input/output rate
grt <- recg$input[recg$ks == min(recg$ks)[1]][1]

# turnover period for debris-flow inputs
tp <- 645.3447740909093
# debris-flow events per turnover period
et <- tp * drt
# debris-flow events remaining
er <- rec$n[rec$ks == min(rec$ks)[1]]

# fluvial events per turnover period
fet <- tp * frt
# fluvial events remaining
fer <- recf$n[recf$ks == min(recf$ks)[1]][1]

# gravel events per turnover period
get <- tp * grt
# gravel events remaining
ger <- recg$n[recg$ks == min(recg$ks)[1]][1]

# total events per turnover
te <- et + fet + get
# total events remaining
tr <- er + fer + ger

# resident fraction
# events remaining divided by events per turnover period
res_fr <- tr / te

# average total valley fill m3
afill <- (vfill + vfillh) / 2

# fill per remaining event
fe <- sfill / tr
# average volume per year flux
af <- fe * te / tp
# denudation rate
dr <- af / (ein * 1e6)

# max denudation rate
mdr <- .000268
# max denudation volume per year flux
mdf <- mdr * ein * 1e6
# time to produce valley storage
mst <- sfill / mdf

# min denudation rate
ndr <- .00007
# min denudation volume per year flux
ndf <- ndr * ein * 1e6
# time to produce valley storage
nst <- sfill / ndf

# from freuh and lancaster
# lower bear creek flux rate
lbf <- sfill / 979 / (ein * 1e6)
# cedar creek fan
ccf <- sfill / 1140 / (ein * 1e6)
# golden ridge creek
grf <- sfill / 1220 / (ein * 1e6)

options(scipen = 6)


rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}

# match site ids between radiocarbon and transect datasets
# note which sites are named differently
crks$sid[!crks$sid %in% charcoal$site_id]
charcoal$site_id[grep('BC', charcoal$site_id)]
crks$sid[grep('BC-18B', crks$sid)] <- 'BC-18'
crks$sid[grep('BC-20B', crks$sid)] <- 'BC-20'
charcoal$site_id[grep('LK', charcoal$site_id)] %>% sort
crks$sid[grep('LK_118R', crks$sid)] <- 'LK_118r'
crks$sid[grep('LK_57b', crks$sid)] <- 'LK_57'
crks$sid[grep('LK_0.42Q', crks$sid)] <- 'LK_0.42q'
crks$sid[grep('LK_140', crks$sid)] <- 'LK140'
crks$sid[grep('LK_36O', crks$sid)] <- 'LK_36o'
crks$sid[grep('LK_0.16L', crks$sid)] <- 'LK_0.16l'
crks$sid[grep('LK_0.16M', crks$sid)] <- 'LK_0.16m'
crks$sid[grep('LK_151N', crks$sid)] <- 'LK_151n'
crks$sid[grep('LK_187L', crks$sid)] <- 'LK_187l'
crks$sid[grep('LK_145F', crks$sid)] <- 'LK_145f'
crks$sid[grep('LK_145M', crks$sid)] <- 'LK_145m'
crks$sid[grep('LK_145O', crks$sid)] <- 'LK_145o'
crks$sid[grep('LK_145P', crks$sid)] <- 'LK_145p'
charcoal$site_id[grep('DFK', charcoal$site_id)] %>% sort
crks$sid[grep('DFK_42', crks$sid)] <- 'DFK_142'
crks$sid[grep('DFK_110B', crks$sid)] <- 'DFK_110b'
crks$sid[grep('DFK_140L', crks$sid)] <- 'DFK_140l'
crks$sid[grep('DFK_140M', crks$sid)] <- 'DFK_140m'
crks$sid[grep('DFK_140N', crks$sid)] <- 'DFK_140n'
crks$sid[grep('DFK_140P', crks$sid)] <- 'DFK_140p'
crks$sid[grep('DFK_140Q', crks$sid)] <- 'DFK_140q'
crks$sid[grep('DFK_140R', crks$sid)] <- 'DFK_140r'
crks$sid[grep('DFK_140S', crks$sid)] <- 'DFK_140s'
crks$sid[grep('DFK_184B', crks$sid)] <- 'DFK_184b'
crks$sid[grep('DFK_280B', crks$sid)] <- 'DFK_280b'

# look up optimized delivery probability and streampower coefficient k for charcoal samples
crks <- creeks_radio
crks_so <- rater1(creeks, .73, .73)
crks$dp <- 0
crks$type <- 0
for (i in 1:nrow(crks)) {
  crks$dp[i] <- crks_so$dp[crks_so$NODE_ID == crks$node_ids[i]]
  crks$k[i] <- crks_so$k[crks_so$NODE_ID == crks$node_ids[i]]
  crks$type[i] <- charcoal$facies[charcoal$site_id == crks$sid[i]]
}





crks$ldp <- log(crks$dp)
md_ldp <- median(crks$ldp)
ad_ldp <- crks$ldp - md_ldp
rng_ldp <- range(ad_ldp)[2]-range(ad_ldp)[1]
wt_dp <- lapply(ad_ldp, function(x) x / rng_ldp) %>% unlist
wt_dp <- -wt_dp

rng_k <- range(crks$k)[2] - range(crks$k)[1]
md_k <- median(crks$k)
wt_k <- lapply(crks$k, function(x) (x - md_k) / rng_k) %>% unlist
wt_k <- -wt_k

df_wt <- sum(wt_dp[crks$type == 'DF'] + 1)
ff_wt <- sum(wt_k[crks$type == 'FF'] + 1)
fg_wt <- sum(wt_k[crks$type == 'FG'] + 1)

# weighted total events remaining
wtr <- df_wt + ff_wt + df_wt

# scale sampled events to total expected events
wt_sc <- tr / nrow(crks) * wtr

# fill per weighted remaining event
wfe <- sfill / wt_sc
# weighted average volume per year flux
waf <- wfe * te / tp
# weighted denudation rate
wdr <- waf / (ein * 1e6)

# KS test
png('ksfit_deposits.png', height = 20, width = 24, units = 'cm', res = 300)
plot(rec$input, rec$ks, pch = 20, col = get_palette('charcoal'),
     xlab = 'input/output rate', ylab = 'Kolmogorov-Smirnov test')
points(recg$input, recg$ks, col = get_palette('ocean'), pch = 20)
points(recf$input, recf$ks, col = get_palette('crimson'), pch = 20)
abline(v = drt, lty = 2, col = get_palette('charcoal', .5) )
abline(v = frt, lty = 2, col = get_palette('crimson', .5) )
abline(v = grt, lty = 2, col = get_palette('ocean', .5) )
legend('topright', legend = c('debris flows', 'fines', 'gravels'),
       fill = get_palette(c('charcoal', 'crimson', 'ocean'), .8))
dev.off()

# print minimum test values
# ad, chi^2, kp, ks debris-flows
c(
  rec$input[rec$ad == min(rec$ad)[1]],
  rec$input[rec$chs == min(rec$chs)[1]],
  rec$input[rec$kp == min(rec$kp)[1]],
  rec$input[rec$ks == min(rec$ks)[1]]
)
# ad, chi^2, kp, ks fluvial fines
c(
  recf$input[recf$ad == min(recf$ad)[1]],
  recf$input[recf$chs == min(recf$chs)[1]],
  recf$input[recf$kp == min(recf$kp)[1]],
  recf$input[recf$ks == min(recf$ks)[1]][1]
)
# ad, chi^2, kp, ks fluvial gravels
c(
  recg$input[recg$ad == min(recg$ad)[1]][1],
  recg$input[recg$chs == min(recg$chs)[1]][1],
  recg$input[recg$kp == min(recg$kp)[1]][1],
  recg$input[recg$ks == min(recg$ks)[1]][1]
)



plot(wt_dp)
plot(wt_k)
plot(crks$k)
crks$wk <- crks$k + abs(min(crks$k))
crks$wk <- crks$wk / max(crks$wk)

plot(crks$wk)
crks_so %>% names

sids <- charcoal$site_id
nrow(charcoal[grep("BC", charcoal$site_id)])
sids <- sids[!sids %in% charcoal$site_id[grep("BC", charcoal$site_id)]]
nrow(charcoal[grep("CC", charcoal$site_id)])
sids <- sids[!sids %in% charcoal$site_id[grep("CC", charcoal$site_id)]]
nrow(charcoal[grep("GRC", charcoal$site_id)])
sids <- sids[!sids %in% charcoal$site_id[grep("GRC", charcoal$site_id)]]
nrow(charcoal[grep("T8", charcoal$site_id)])
sids <- sids[!sids %in% charcoal$site_id[grep("T8", charcoal$site_id)]]
nrow(charcoal[grep("DFK", charcoal$site_id)])
sids <- sids[!sids %in% charcoal$site_id[grep("DFK", charcoal$site_id)]]
nrow(charcoal[grep("UK", charcoal$site_id)])
sids <- sids[!sids %in% charcoal$site_id[grep("UK", charcoal$site_id)]]
nrow(charcoal[grep("LK", charcoal$site_id)])
sids <- sids[!sids %in% charcoal$site_id[grep("LK", charcoal$site_id)]]

# odds of flux selection
pf <- te / (te + tr)

fish <- function(rt, t, k) {
  (rt*t)^k*exp(-rt*t) / factorial(k)
}

df_cdf <- emp_cdf(charcoal$mn[charcoal$facies == 'DF'])
df_val <- emp_cdf(charcoal$mn[charcoal$facies == 'DF'])[,1]
plot(emp_cdf(charcoal$mn[charcoal$facies == 'DF']))
plot(pex)

t <- seq(0, 37, length.out = 4430)
# probability of exiting single reservoir within turnover period t (k is exit)
pex_pmf <- (fish(1-pf, t, 1)/sum(fish(1-pf,t,1)))
index <- char_pmfs %>% rownames %>% as.numeric %>% rev
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
tb <- seq(0, 1, length.out = 4430)
pex_b <- (fish(pf, t, 1)/sum(fish(pf,t,1)))
plot(tb * tp, pex_b)
plot(t * tp, pex_pmf)
lines(tb * tp, pex_b)

index <- char_pmfs %>% rownames %>% as.numeric %>% rev
cfl <- rep(0, length(index))
cfl[1] <- sum(pex_pmf[t * tp < index[1]])
for (i in 2:length(index)) {
  cfl[i] <- sum(pex_pmf[(t * tp) > index[i-1] & (t * tp) <= index[i]])
}
cfl <- convo_add(pex_pmf, df_pmf, index)
cflb <- rep(0, length(index))
cflb[1] <- sum(pex_b[tb * tp < index[1]])
for (i in 2:length(index)) {
  cflb[i] <- sum(pex_b[(tb * tp) > index[i-1] & (tb * tp) <= index[i]])
}

cflb <- convo(df_pmf, cflb, index)
tfl <- cfl * (tr/(te+tr)) + cflb * (te/(te+tr))
plot(index, cumsum((cflb + cfl)/2))
plot(index, cumsum(df_pmf))
plot(index, cumsum(pex_b))
lines(index, cumsum(df_pmf))

# transform cdfs using fft
tdf <- fft(df_cdf) / length(df_cdf)
tpex <- fft(pex) / length(pex)
# add tranformed fluvial inheritance to debris flows
tfl <- tdf + tpex
# invert transform
ifl <- fft(tfl, inverse = T)
plot(ifl)



plot(tdf)
k <- 1:10
png('exp_fluvial_inher.png', height = 20, width = 24, units = 'cm', res = 300)
plot(t*tp, cumsum(pex_pmf),
     xlab = 'expected fluvial inheritence [years]',
     ylab = 'CDF', xlim = c(0, 4000),
     type = 'l', lwd = 2, col = get_palette('ocean'))
points(t*tp, cumsum(pex_pmf), pch = 20, col = get_palette('charcoal', .01))
lines(emp_cdf(charcoal$mn[charcoal$facies == 'FF']))
lines(emp_cdf(charcoal$mn[charcoal$facies == 'FG']))
dev.off()

plot(t, pex)
for (i in k) {
  lines(t, fish(pf, t, i))
}
cumsum(fish(pf, t, 2))/sum(fish(pf,t,2))

factorial(1000)

event_stats <- function(drt, frt, grt, der, fer, ger, tp, fill, so) {
  # debris-flow events per turnover period
  det <- tp * drt
  # fluvial events per turnover period
  fet <- tp * frt
  # gravel events per turnover period
  get <- tp * grt
  # total events per turnover
  te <- det + fet + get
  # total events remaining
  tr <- der + fer + ger
  # fill per remaining event
  fe <- fill / tr
  # average volume per year flux
  af <- fe * te / tp
  # denudation rate
  dr <- af / (ein * 1e6)

  # create debris-flow weights
  # log delivery probability p
  ldp <- log(so$dp)
  # median log p
  md_ldp <- median(ldp)
  # distance from median p in log space
  ad_ldp <- ldp - md_ldp
  # range of variation in log space
  rng_ldp <- range(ad_ldp)[2]-range(ad_ldp)[1]
  # relative distance, scaled by range
  wt_dp <- lapply(ad_ldp, function(x) x / rng_ldp) %>% unlist
  # negated to indicate increasing significance as p goes down
  wt_dp <- -wt_dp

  # create fluvial weights
  rng_k <- range(so$k)[2] - range(so$k)[1]
  md_k <- median(so$k)
  wt_k <- lapply(so$k, function(x) (x - md_k) / rng_k) %>% unlist
  wt_k <- -wt_k

  df_wt <- sum(wt_dp[so$type == 'DF'] + 1)
  ff_wt <- sum(wt_k[so$type == 'FF'] + 1)
  fg_wt <- sum(wt_k[so$type == 'FG'] + 1)

  # weighted total events remaining
  wtr <- df_wt + ff_wt + fg_wt
  # scale sampled events to total expected events
  wt_sc <- tr / nrow(so) * wtr
  # fill per weighted remaining event
  wfe <- sfill / wt_sc
  # weighted average volume per year flux
  waf <- wfe * te / tp
  # weighted denudation rate
  wdr <- waf / (ein * 1e6)
  df <- data.frame(fe = fe,
                   af = af,
                   dr = dr,
                   wfe = wfe,
                   waf = waf,
                   wdr = wdr)
  df
}

# anderson-darling
# debris-flow input/output rate
drt <- rec$input[rec$ad == min(rec$ad)[1]]
# fines input/output rate
frt <- recf$input[recf$ad == min(recf$ad)[1]][1]
# gravels input/output rate
grt <- recg$input[recg$ad == min(recg$ad)[1]][1]

# debris-flow events remaining
der <- rec$n[rec$ad == min(rec$ad)[1]][1]
# fines events remaining
fer <- recf$n[recf$ad == min(recf$ad)[1]][1]
# gravel events remaining
ger <- recg$n[recg$ad == min(recg$ad)[1]][1]

# turnover period for debris-flow inputs
tp <- 318

event_stats(drt, frt, grt, der, fer, ger, tp, sfill, crks)

# chi-squared
# debris-flow input/output rate
drt <- rec$input[rec$chs == min(rec$chs)[1]]
# fines input/output rate
frt <- recf$input[recf$chs == min(recf$chs)[1]][1]
# gravels input/output rate
grt <- recg$input[recg$chs == min(recg$chs)[1]][1]

# debris-flow events remaining
der <- rec$n[rec$chs == min(rec$chs)[1]][1]
# fines events remaining
fer <- recf$n[recf$chs == min(recf$chs)[1]][1]
# gravel events remaining
ger <- recg$n[recg$chs == min(recg$chs)[1]][1]

# turnover period for debris-flow inputs
tp <- 293

event_stats(drt, frt, grt, der, fer, ger, tp, sfill, crks)

# kuiper
# debris-flow input/output rate
drt <- rec$input[rec$kp == min(rec$kp)[1]]
# fines input/output rate
frt <- recf$input[recf$kp == min(recf$kp)[1]][1]
# gravels input/output rate
grt <- recg$input[recg$kp == min(recg$kp)[1]][1]

# debris-flow events remaining
der <- rec$n[rec$kp == min(rec$kp)[1]][1]
# fines events remaining
fer <- recf$n[recf$kp == min(recf$kp)[1]][1]
# gravel events remaining
ger <- recg$n[recg$kp == min(recg$kp)[1]][1]

# turnover period for debris-flow inputs
tp <- 191

event_stats(drt, frt, grt, der, fer, ger, tp, sfill, crks)

# kolmogorov-smirnov
# debris-flow input/output rate
# Debris-flow transit time quantiles are
# [64.48734880227971, 284.3732106558563, 598.8059856023998, 1348.0931486710629, 6899.373207597647]
drt <- rec$input[rec$ks == min(rec$ks)[1]]
# fines input/output rate
frt <- recf$input[recf$ks == min(recf$ks)[1]][1]
# gravels input/output rate
grt <- recg$input[recg$ks == min(recg$ks)[1]][1]

# debris-flow events remaining
der <- rec$n[rec$ks == min(rec$ks)[1]][1]
# fines events remaining
fer <- recf$n[recf$ks == min(recf$ks)[1]][1]
# gravel events remaining
ger <- recg$n[recg$ks == min(recg$ks)[1]][1]

# turnover period for debris-flow inputs
tp <- 598.8059856023998

ks_stats <- event_stats(drt, frt, grt, der, fer, ger, tp, sfill, crks)

df_stereo <- data.table::fread('/home/erik/output/df_stereo.csv') %>% unlist
df_stereo_ad <- data.table::fread('/home/erik/output/df_stereo_ad.csv') %>% unlist
df_stereo_ch <- data.table::fread('/home/erik/output/df_stereo_ch.csv') %>% unlist
df_stereo_kp <- data.table::fread('/home/erik/output/df_stereo_kp.csv') %>% unlist
# debris-flow transit times at selected fits
df_transit_ad <- data.table::fread('/home/erik/output/df_transits_ad.csv') %>% unlist
df_transit_ch <- data.table::fread('/home/erik/output/df_transits_ch.csv') %>% unlist
df_transit_kp <- data.table::fread('/home/erik/output/df_transits_kp.csv') %>% unlist
df_transit_ks <- data.table::fread('/home/erik/output/df_transits_ks.csv') %>% unlist
ff_stereo <- data.table::fread('/home/erik/output/ff_stereo.csv') %>% unlist
fg_stereo <- data.table::fread('/home/erik/output/fg_stereo.csv') %>% unlist

# subset debris-flow charcoal means in bear and knowles
charsub <- charcoal$mn[charcoal$site_id %in% crks$sid[crks$creek_name %in% c('knowles', 'bear')] &
                         charcoal$facies == 'DF']

png('df_stero_fit.png', height = 20, width = 24, units = 'cm', res = 300)
plot(emp_cdf(charsub),
     xlab = 'charcoal age [y]', ylab = 'CDF',
     type = 'l', lwd = 2, col = get_palette('charcoal', .01))
lines(emp_cdf(df_stereo), lwd = 2, col = get_palette('crimson'))
points(emp_cdf(df_stereo), pch = 20, col = get_palette('crimson', .2))
lines(emp_cdf(df_stereo_ad), lwd = 2, col = get_palette('ocean'))
points(emp_cdf(df_stereo_ad), pch = 20, col = get_palette('ocean', .2))
lines(emp_cdf(df_stereo_ch), lwd = 2, col = get_palette('gold'))
points(emp_cdf(df_stereo_ch), pch = 20, col = get_palette('gold', .2))
lines(emp_cdf(df_stereo_kp), lwd = 2, col = get_palette('violet'))
points(emp_cdf(df_stereo_kp), pch = 20, col = get_palette('violet', .2))
lines(emp_cdf(charsub), col = get_palette('charcoal'))
points(emp_cdf(charsub),
       pch = 20, col = get_palette('charcoal'))
legend('bottomright', legend = c('observed', 'anderson-darling',
                                 'chi-squared', 'kuiper', 'kolmogorov-smirnov'),
       fill = get_palette(c('charcoal', 'ocean', 'gold', 'violet', 'crimson'), .8))
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

png('df_stero_fit.png', height = 20, width = 24, units = 'cm', res = 300)
plot(emp_cdf(charcoal$mn[charcoal$facies == 'DF']),
     xlab = 'charcoal age', ylab = 'CDF',
     type = 'l', lwd = 2, col = get_palette('charcoal'))
lines(emp_cdf(df_stereo), lwd = 2, col = get_palette('crimson'))
points(emp_cdf(charcoal$mn[charcoal$facies == 'DF']),
       pch = 20, col = get_palette('charcoal'))
points(emp_cdf(df_stereo), pch = 20, col = get_palette('crimson', .2))
legend('bottomright', legend = c('observed', 'fit'),
       fill = get_palette(c('charcoal', 'crimson'), .8))
dev.off()

png('fl_stereo_fit.png', height = 20, width = 24, units = 'cm', res = 300)
plot(emp_cdf(charcoal$mn[charcoal$facies == 'FG']),
     xlab = 'charcoal age', ylab = 'CDF',
     type = 'l', lwd = 2, col = get_palette('ocean'))
lines(emp_cdf(fg_stereo), lwd = 2, col = get_palette('sky'))
points(emp_cdf(charcoal$mn[charcoal$facies == 'FG']),
       pch = 20, col = get_palette('ocean'))
points(emp_cdf(fg_stereo), pch = 20, col = get_palette('sky', .2))
lines(emp_cdf(ff_stereo), lwd = 2, col = get_palette('rose'))
points(emp_cdf(ff_stereo), pch = 20, col = get_palette('rose', .2))
lines(emp_cdf(charcoal$mn[charcoal$facies == 'FF']),
      lwd = 2, col = get_palette('crimson'))
points(emp_cdf(charcoal$mn[charcoal$facies == 'FF']),
       pch = 20, col = get_palette('crimson'))
legend('bottomright', legend = c('obs gravels', 'fit gravels', 'obs fines', 'fit fines'),
       fill = get_palette(c('ocean', 'sky', 'crimson', 'rose'), .8))
dev.off()

