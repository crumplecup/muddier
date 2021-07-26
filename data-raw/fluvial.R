library(data.table)
library(magrittr)
library(parallel)
library(plot3D)
library(sp)
library(muddier)
options(mc.cores=detectCores())
set.seed(10101)



# subset sites from bear and knowles
# remove older samples from sites with multiples

site <- creeks_radio$node_ids
counts <- 0
for (i in seq_along(site)) {
  counts[i] <- length(site[site == site[i]])
}
mults <- site[counts > 1] %>% factor %>% levels

crks <- creeks_radio[!creeks_radio$node_ids %in% mults, ]

for (i in seq_along(mults)) {
  sub <- creeks_radio[creeks_radio$node_ids == mults[i], ]
  crks <- rbind(crks, sub[1,])
}

crks <- crks[crks$creek_name != 'cedar', ]
rm(sub, counts, i, mults, site)




# sort charcoal pmfs by facies and remove oversamples
sid <- char_pmfs %>% colnames
pmfs <- t(char_pmfs)

# unique deposit sites
ct <- charcoal[ , .N, by = c('facies', 'family')]

# minimum sample age from each deposit
ids <- 0
for (i in 1:nrow(ct)) {
  pool <- charcoal[charcoal$facies == ct$facies[i] & charcoal$family == ct$family[i]]
  ids[i] <- pool$site_id[pool$mn == min(pool$mn)]
}

# pmf of each site by id
ct$ids <- ids
ar <- array(0, c(length(ids), ncol(pmfs)))
for (i in seq_along(ids)) {
  ar[i, ] <- pmfs[sid == ids[i], ]
}

index <- char_pmfs %>% rownames %>% as.numeric

rownames(ar) <- ids
colnames(ar) <- index

# subset by facies
df <- ar[ct$facies == 'DF', ]
ff <- ar[ct$facies == 'FF', ]
fg <- ar[ct$facies == 'FG', ]

nrow(df)
nrow(ff)
nrow(fg)


# interarrival times for debris flows

# order by weighted mean
index <- char_pmfs %>% rownames %>% as.numeric %>% sort
dfmn <- apply(df, 1, function (x) weighted.mean(rev(index), x))
itmn <- 0
for (i in 2:length(dfmn)) {
  itmn[i] <- dfmn[i] - dfmn[i-1]
}

ffmn <- apply(ff, 1, function (x) weighted.mean(rev(index), x))
ffitmn <- 0
for (i in 2:length(ffmn)) {
  ffitmn[i] <- ffmn[i] - ffmn[i-1]
}

fgmn <- apply(fg, 1, function (x) weighted.mean(rev(index), x))
fgitmn <- 0
for (i in 2:length(fgmn)) {
  fgitmn[i] <- fgmn[i] - fgmn[i-1]
}


# remove GRC and Cedar creek samples for now

dfids <- rownames(df)
ffids <- rownames(ff)
fgids <- rownames(fg)

rems <- c(grep('GRC', rownames(df)), grep('CC', rownames(df)))
dfs <- df[-rems,]
dfids <- dfids[-rems]
rems <- c(grep('GRC', rownames(ff)), grep('CC', rownames(ff)), grep('T8', rownames(ff)))
ffs <- ff[-rems,]
ffids <- ffids[-rems]
rems <- c(grep('GRC', rownames(fg)), grep('CC', rownames(fg)), grep('T8', rownames(fg)))
fgs <- fg[-rems,]
fgids <- fgids[-rems]
rm(rems, ids)


pal <- get_palette(c('ocean', 'forest', 'gold', 'crimson'))
plot(emp_cdf(fgmn), pch = 20, col = pal[1])
points(emp_cdf(ffmn), pch = 20, col = pal[3])
points(emp_cdf(dfmn), pch = 20, col = pal[2])

png('facies_interarrival.png', height = 17, width = 21, units = 'cm', res = 300)
pal <- get_palette(c('crimson', 'gold', 'ocean', 'charcoal'))
plot(fgmn, fgitmn, log = 'x', col = pal[3], pch = 20,
     xlim = c(min(dfmn), max(fgmn)),
     xlab = 'charcoal age [years]',
     ylab = 'interarrival time [years]')
points(dfmn, itmn, pch = 20, col = pal[4])
points(ffmn, ffitmn, pch = 20, col = pal[1])
legend('topleft', legend = c('debris flows', 'gravels', 'fines'),
       fill = get_palette(c('charcoal', 'ocean', 'crimson'), .9))
dev.off()



library(lubridate)

flow <- data.table::fread('/media/erik/catacomb/research/siuslaw_flow.csv')
flow$dt <- mdy_hm(flow$datetime)
flow$yr <- year(flow$dt)
flow$mt <- month(flow$dt)
flow$dr <- 0
flow$dr[-1] <- flow$dt[-1] - flow$dt[-length(flow$dt)]
flow$dr[1] <- flow$dr[2]
flow$dr[flow$dr > 1800] <- 1800
flow$dq <- 0
flow$dq[-1] <- flow$discharge_cfs[-1] - flow$discharge_cfs[-length(flow$discharge_cfs)]

is_peak <- function(vec) {
  vec[1] > 0 & vec[2] < 0
}

peaker <- function(vec) {
  ln <- length(vec)
  flag <- vector(length = ln, mode = 'numeric')
  for(i in seq_along(vec)) {
    if (i < ln) {
      flag[i] <- vec[i] > 0 & vec[i+1] < 0
    }
  }
  flag
}

peaks <- peaker(dec20$dq)
peaks
plot(dec20$discharge_cfs[peaks == 1])
length(sum(peaks)) / length(peaks)

library(hydrostats)
library(magrittr)
# subset flows after 2000, due to a gap in coverage, and different gauge times
df <- data.frame(Date = flow$dt[flow$yr > 2000] %>% as.POSIXct,
                 Q = flow$discharge_cfs[flow$yr > 2000],
                 dt = flow$dt[flow$yr > 2000])
df <- df[!is.na(df$Q), ]
head(df)
ann.cv(df)
baseflows(df)
baseflows(df, ts = 'annual')
Colwells(df)
CTF(df)
daily.cv(df)
flood.length.max(df, 40000)
high.spell.lengths(df)
high.spells(df, threshold = 30000)
ps <- partial.series(df, ari = 1, series = T, duration = T, plot = T)
seasonality(df, monthly = T)

# subset flows above a 'critical' minimum threshold
crit_flow <- 18400
hiflo <- data.frame(
  Date = flow$dt[flow$yr > 2000 & flow$discharge_cfs >= crit_flow] %>% as.POSIXct,
  Q = flow$discharge_cfs[flow$yr > 2000 & flow$discharge_cfs >= crit_flow])
hi_cdf <- emp_cdf(hiflo$Q)
hi_df <- hi_cdf %>% as.data.frame
hi_mod <- lm(y ~ log(x), data = hi_df)
hi_pred <- predict(hi_mod, newdata = hi_df)
summary(hi_mod)
plot(hi_cdf, log = 'x')
lines(hi_df$x, hi_pred/max(hi_pred))

ps <- partial.series(df[flow$yr > 2000, ], series = T)
psr <- ps$p.series %>% as.data.table
psr <- psr[order(Date), ]
psr
psd <- psr$Date %>% as_datetime
psx <- psd[-1] - psd[-length(psd)]
plot(psr$Q[-1], sort(psx))

length(psx) / sum(psx) %>% as.numeric * 365

dif <- sum(psx)

sius_ca <- 588 / 0.38610 # sq miles from USGS gauge site converted to km2
sa_flow <- flow[!is.na(flow$discharge_cfs), ]
sa_flow <- sa_flow[sa_flow$yr > 2000, ]
sa_flow$Q <- sa_flow$discharge_cfs * (ein / sius_ca)
sa_flow$rq <- sa_flow$Q / sum(sa_flow$Q)
sa_flow$q <- sa_flow$Q * 0.028316847

# look up optimized delivery probability and streampower coefficient k for charcoal samples
crks <- creeks_radio
crks_so <- rater1(creeks, .73, .73)
# convert lengths from to-mouth to hipchain from study area
crks_so$hip <- crks_so$ToMouth_km
crks_so$hip[crks_so$creek_name == 'knowles'] <- crks_so$hip[crks_so$creek_name == 'knowles'] - min(crks_so$hip[crks_so$creek_name == 'knowles'])

# subset bear and knowles
bk <- crks_so[crks_so$creek_name %in% c('bear', 'knowles'), ]
hip_cdf <- emp_cdf(bk$hip)

map_cdf <- function(obs) {
  cdf <- emp_cdf(obs)
  unlist(lapply(obs, function(x) max(cdf[cdf[ ,1] <= x , 2])))
}

plot(density(bk$hip))
nrow(bk)
dim(emp_cdf(bk$hip))
bk$cdf <- map_cdf(bk$hip)
plot(to_pmf(emp_cdf(bk$hip)[,2]))
plot(bk$hip, hip_cdf1)

# subset bear then knowles
br <- bk[bk$creek_name == 'bear', ]
kn <- bk[bk$creek_name == 'knowles', ]

# from km to m
# switch orientation from mouth to initiation point
br$hp <- (max(br$hip) - br$hip) * 1000
kn$hp <- (max(kn$hip) - kn$hip) * 1000
library(magrittr)
# cumulative lengths
br$cm <- br$hp %>% rev %>% cumsum %>% rev
kn$cm <- kn$hp %>% rev %>% cumsum %>% rev
# upstream count
br$ct <- nrow(br):1
kn$ct <- nrow(kn):1
# average length traveled per segment
br$alt <- br$cm / br$ct
kn$alt <- kn$cm / kn$ct
# fluvial inheritance [y m^-1]
mn_d <- mean(c(br$hip, kn$hip)) # mean travel distance [km]
br$fi <- br$alt / (mn_d * 1000 / tp)
kn$fi <- kn$alt / (mn_d * 1000 / tp)

plot(br$hip*1000)
lines(br$hp)
plot(kn$hp, kn$alt, pch = 20, col = get_palette('ocean'))
points(kn$hp[kn$NODE_ID %in% crks$node_ids], kn$alt[kn$NODE_ID %in% crks$node_ids],
       pch = 20, col = get_palette('crimson'))
points(br$hp, br$alt, pch = 20, col = get_palette('ocean'))
points(br$hp[br$NODE_ID %in% crks$node_ids], br$alt[br$NODE_ID %in% crks$node_ids],
       pch = 20, col = get_palette('crimson'))


ch_alt <- c(br$alt[br$NODE_ID %in% crks$node_ids], kn$alt[kn$NODE_ID %in% crks$node_ids])

plot(ch_alt / mn_d, ylab = 'fluvial inheritance (y)')
abline(h = mean(ch_alt / mn_d))
# mean reach length [m]
mn_ln <- mean(bk$hip) * 1000
# mean sediment distance per year [m/y]
mn_d <- mn_ln / tp
# average time from mouth [y]
crks_so$tm <- crks_so$hip * 1000 * mn_d


k <- 1:10
t <- 1:50
r <- 0.5
plot(t, ((r*t)^1*(exp(-r*t)))/factorial(1))
for (i in k) {
  lines(t, ((r*t)^i*(exp(-r*t)))/factorial(i))

}

png('fl_inher.png', height = 17, width = 25, units = 'cm', res = 300)
plot(kn$hp, kn$fi, pch = 20, col = get_palette('ocean'), type = 'l', lwd = 2,
     xlab = 'distance from initiation point [m]',
     ylab = 'average traversal time [y]')
points(kn$hp[kn$NODE_ID %in% crks$node_ids], kn$fi[kn$NODE_ID %in% crks$node_ids],
       pch = 20, col = get_palette('crimson'), cex = .8)
lines(br$hp, br$fi, pch = 20, col = get_palette('sky'), lwd = 2)
points(br$hp[br$NODE_ID %in% crks$node_ids], br$fi[br$NODE_ID %in% crks$node_ids],
       pch = 20, col = get_palette('crimson'), cex = .8)
legend('bottomright', legend = c('bear', 'knowles', 'charcoal'),
       fill = get_palette(c('sky', 'ocean', 'crimson'), .8))
dev.off()

(max(br$fi) + max(kn$fi))/2

dual_cdf <- function(x, y) {
  xl <- length(x)
  yl <- length(y)
  # join the two obs, sort and deduplicate
  z <- sort(unique(c(x, y)))
  x_f <- 0
  y_f <- 0
  for (i in 1:length(z)) {
    x_f[i] <- length(x[x <= z[i]]) / xl
    y_f[i] <- length(y[y <= z[i]]) / yl
  }
  list(x, y)
}

# fit interarrival times to a lambda
iat <- data.table::fread('/home/erik/data/iat.csv')
iat <- iat[order(iat$ia), ]
iat$iat <- 0
for (i in 2:nrow(iat)) {
  iat$iat[i] <- iat$ia[i] - iat$ia[i-1]
}

plot(iat$ia %>% sort)
plot(iat$iat, iat$ia)

xy <- coordinates(crks_so)
crks_so$x <- xy[,1]
crks_so$y <- xy[,2]

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

crks$dp <- 0
crks$type <- 0
crks$hip <- 0
crks$mn <- 0
crks$tm <- 0
for (i in 1:nrow(crks)) {
  crks$dp[i] <- crks_so$dp[crks_so$NODE_ID == crks$node_ids[i]]
  crks$k[i] <- crks_so$k[crks_so$NODE_ID == crks$node_ids[i]]
  crks$hip[i] <- crks_so$hip[crks_so$NODE_ID == crks$node_ids[i]]
  crks$tm[i] <- crks_so$tm[crks_so$NODE_ID == crks$node_ids[i]]
  crks$type[i] <- charcoal$facies[charcoal$site_id == crks$sid[i]]
  crks$mn[i] <- charcoal$mn[charcoal$site_id == crks$sid[i]]
}


plot(crks$tm[crks$type != 'DF'], crks$mn[crks$type != 'DF'] / crks$tm[crks$type != 'DF'],
     log = 'y')
points(crks$tm[crks$type == 'DF'], crks$mn[crks$type == 'DF'] / crks$tm[crks$type == 'DF'],
       pch = 20, col = get_palette('crimson'))

plot(emp_cdf(crks_so$tm[crks$creek_name %in% c('knowles', 'bear')]))
abline(v = mn_ln*1000)

plot(crks$hip[crks$type != 'DF'], crks$mn[crks$type != 'DF'])
abline(lm(mn ~ hip, data = crks))
abline(h = mn_ln)
points(creeks$hip[creeks$NODE_ID %in% creeks_radio$node_ids],
       pch = 20, col = get_palette('crimson'))


# total sediment storage in the study creeks
sedvol <- data.table::fread('/home/erik/data/sedvol.csv')
sfill <- sum(sedvol$`Volume increment (m3)`[sedvol$Reach == 'Cedar']) +
  sum(sedvol$`Volume increment (m3)`[sedvol$Reach == 'Hoffman']) +
  sum(sedvol$`Volume increment (m3)`[sedvol$Reach %in% c('LB', 'UB')]) +
  sum(sedvol$`Volume increment (m3)`[sedvol$Reach %in% c('DFK', 'LK', 'UK')])
# contributing area for study creeks
ein <- max(creeks$contr_area[creeks$creek_name == 'bear']) +
  max(creeks$contr_area[creeks$creek_name == 'cedar']) +
  max(creeks$contr_area[creeks$creek_name == 'hoffman']) +
  max(creeks$contr_area[creeks$creek_name == 'knowles'])
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

# sq miles from USGS gauge site converted to km2
sius_ca <- 588 / 0.38610
# flow scaled to study area
df$Qsa <- df$Q * (ein / sius_ca)
# relative flow
df$rq <- df$Qsa / sum(df$Qsa)
# convert flow from cfs to m3/s
df$q <- df$Qsa * 0.028316847
# duration of flow measure
gage_t <- interval(min(df$dt), max(df$dt)) / years(1)
# q2r for total discharge
# discharge squared
df$q2 <- df$q^2
# relative squared discharge
df$q2r <- df$q2 / sum(df$q2)
# years duration * weighted average flux * relative force applied
df$v <- gage_t * ks_stats$waf * df$q2r
# duration of gauge measure is every 15m
# convert from m3/15m to m3/s
df$vs <- df$v / (60 * 15)
# sediment density of 2650 kg/m3
# convert from m3/s to mg/l
df$mgl <- 2650 * df$vs / df$q * 1000
plot(df$q, df$mgl, log = 'y')


# critical motion thresholds for sand and gravel from wilcock and kenworthy 2002
# sand tau_ri
t_rs <- 7.5
# gravel tau_ri
t_rg <- 12

# 0.84 * D_90 for Oak Creek [mm]
# surface D_90 for Oak Creek is 98 mm
# from wilcock & kenworthy (2002)
# To provide consistency with the laboratory data, we
# set k_s = 0.84 * D_90, a value back calculated from the sidewall-
# corrected values of t in the flume experiments
d90 <- 0.84 * 98
# diameter of average gravel [mm] (surface)
dg <- 53
# diameter of average sand [mm] (subsurface)
ds <- 1.2

# from wilcock etal 2009
# shear stress due to bed grains only
# using the Manning-Strickler formula
# tau prime = rho * g * (0.013)^(3/2) * (S * D)^(1/4) * U^(3/2)
tau1 <- 1000 * 9.81 * (0.013)^(3/2) * (0.074 * d90)^(1/4) * df$Qkm2^(3/2)
png('bed_shear_stress.png', height = 17, width = 21, units = 'cm', res = 300)
plot(df$dt, tau1, type = 'l', lwd = 2, col = get_palette('ocean'),
     xlab = 'year', ylab = 'bed shear stress (Pa)', log = 'y')
dev.off()

# dimensionless conversion for tau_star
to_tau_star <- function(tau, d) {
  tau / (1650 * 9.81 * d)
}

# dimensionless conversion for q_star
to_qstar <- function(q, fr, d) {
  q / (fr * sqrt(1.650 * 9.81 * d^3))
}

# annual sediment flux per km2 in m3
fkm2 <- ks_stats$waf / ein



# flow per km2
df$Qkm2 <- df$q / ein
# k38p2026 is 2m wide with 0.4m vertical banks and \approx 1km contr area
# flow height per meter of width
df$h <- df$Qkm2 / 2
# wetted perimeter [m]
df$P <- df$h * 2 + 2
# hydraulic radius [m]
df$R <- df$Qkm2 / df$P
# shear stress at p2026
# [kg/m3] * [m/s^2] * [m] = kg m^-1 s^-2
# in pascals (drop)
# shear stress = density of water * gravity * hydraulic radius * slope
df$t0 <- 1000 * 9.81 * df$R * 0.074
plot(df$dt, df$t0, type = 'l', lwd = 2, col = get_palette('ocean'),
     xlab = 'year', ylab = 'shear stress (Pa)')




png('flux_q.png', height = 23, width = 30, units = 'cm', res = 300)
plot(sa_flow$dt, sa_flow$q, type = 'l', log = 'y',
     ylim = c(min(sa_flow$v / (60 * 15)), max(sa_flow$q)),
     xlab = 'year', ylab = 'Flux [m3/s]',
     col = get_palette('ocean', .5))
lines(sa_flow$dt, sa_flow$v / (60 * 15), col = get_palette('crimson', .5))
legend('left', legend = c('Total Discharge', 'Sediment'),
       fill = get_palette(c('ocean', 'crimson'), .8))
dev.off()

gage_t <- interval(min(sa_flow$dt), max(sa_flow$dt)) / years(1)
sa_flow$v <- gage_t * ks_stats$waf * sa_flow$rq

png('vol_by_q.png', height = 23, width = 30, units = 'cm', res = 300)
plot(sa_flow$q, sa_flow$v / (60 * 15), type = 'l', lwd = 2,
     xlab = 'Total Discharge [m3/s]', ylab = 'Sediment Flux [m3/s]',
     col = get_palette('charcoal', .8))
points(sa_flow$q, sa_flow$v / (60 * 15), pch = 20,
       col = get_palette('coral', .025))
legend('bottomright', legend = c('observed', 'fit'),
       fill = get_palette(c('coral', 'charcoal'), .8))
dev.off()

png('flux_by_q.png', height = 23, width = 30, units = 'cm', res = 300)
plot(sa_flow$dt, sa_flow$v / (60 * 15), type = 'l', col = get_palette('ocean'),
     xlab = 'date', ylab = 'sediment flux [m3/s]')
dev.off()
png('flux_bedload.png', height = 23, width = 30, units = 'cm', res = 300)
plot(sa_flow$dt, (sa_flow$v / (60 * 15)) / sa_flow$q, type = 'l', col = get_palette('ocean'),
     xlab = 'date', ylab = 'sediment flux [m3/s]')
dev.off()
plot(emp_cdf(sa_flow$v), type = 'l', col = get_palette('ocean'))
ks_stats$waf


q_ord <- sa_flow$q %>% sort
q_cum <- cumsum(q_ord)
plot(q_cum)
plot(q_ord, q_cum)

plot(sa_flow$dt, sa_flow$Q)
plot(sa_flow$rq, sa_flow$Q)

# subset flows above the indicated recurrence interval
df <- data.frame(Date = flow$dt[flow$yr > 2000] %>% as.POSIXct,
                 Q = flow$discharge_cfs[flow$yr > 2000])
thr <- partial.series(df, series = T, ari = 0.666)
thresh <- thr$flow.threshold

# map sediment proportionally to subset discharge
sa_sub <- sa_flow[sa_flow$discharge_cfs >= thresh, ]
sa_sub$rq <- sa_sub$Q / sum(sa_sub$Q)
sa_sub$v <- gage_t * ks_stats$waf * sa_sub$rq
# discharge squared
sa_sub$q2 <- sa_sub$q^2
# relative squared discharge
sa_sub$q2r <- sa_sub$q2 / sum(sa_sub$q2)
# sediment density of 2650 kg/m3
sa_sub$v <- gage_t * ks_stats$waf * sa_sub$q2r
sa_sub$mgl <- 2650 * sa_sub$v / (60 * 15) / sa_sub$q * 1000
plot(sa_sub$q, sa_sub$mgl, log = 'y')


sa_flow <- sa_flow[order(sa_flow$q), ]

png('vol_by_events.png', height = 23, width = 30, units = 'cm', res = 300)
plot(sa_flow$q, sa_flow$mgl, type = 'l', lwd = 2,
     xlab = 'Total Discharge [m3/s]', ylab = 'Sediment Concentration [mg/l]',
     col = get_palette('charcoal', .8),
     ylim = c(min(sa_flow$mgl), max(sa_sub$mgl)), log = 'y')
points(sa_flow$q, sa_flow$mgl, pch = 20,
       col = get_palette('sky', .025))
lines(sa_sub$q, sa_sub$mgl, lwd = 2,
      col = get_palette('charcoal'))
points(sa_sub$q, sa_sub$mgl, pch = 20,
      col = get_palette('coral', .05))
legend('topleft', legend = c('fit to discharge', 'fit to events'),
       fill = get_palette(c('sky', 'coral'), .8))
dev.off()

png('flux_q.png', height = 23, width = 30, units = 'cm', res = 300)
plot(sa_flow$dt, sa_flow$q, type = 'l', log = 'y',
     xlab = 'year', ylab = 'Flux [m3/s]',
     col = get_palette('ocean', .5),
     ylim = c(min(sa_flow$vs), max(sa_flow$q)))
lines(sa_flow$dt, sa_flow$vs, col = get_palette('crimson', .5))
legend('bottomleft', legend = c('Total Discharge', 'Sediment'),
       fill = get_palette(c('ocean', 'crimson'), .8))
dev.off()



qdf <- data.frame(cdf = q_cdf[,2], q = q_cdf[,1])
q_mod <- lm('cdf ~ log(q)', data = qdf[qdf$q > crit_flow, ])
summary(q_mod)
q_pred <- predict(q_mod, newdata = qdf[qdf$q > crit_flow, ])
plot(qdf$q[qdf$q > crit_flow], qdf$cdf[qdf$q > crit_flow], log = 'x')
lines(qdf$q[qdf$q > crit_flow], qdf$cdf[qdf$q > crit_flow],
      col = get_palette('crimson'))
lines(qdf$q[qdf$q > crit_flow], q_pred/max(q_pred))
plot(q_pred/max(q_pred))

lines(q_pred, col = get_palette('crimson'))
u <- 1
z <- .07
y <- ( u * (1:50000)^z)
plot(1:50000, y)
lines(qdf$q, qdf$cdf)

length(flow$discharge_cfs[flow$discharge_cfs > q_cdf[1000,1]])
qdf$n <- lapply(qdf$q, function(x) length(flow$discharge_cfs[flow$discharge_cfs >= x])) %>% unlist
plot(qdf$q, qdf$n)
max_d <- max(qdf$q)
qdf$rq <- qdf$q / max_d
plot(qdf$rq)
plot(qdf$n, qdf$rq, log = 'x')
N <- nrow(flow)
q_mod <- lm('log(n) ~ log(rq)', data = qdf)
summary(q_mod)
q_pred <- predict(q_mod, newdata = qdf) %>% exp
plot(qdf$q, qdf$n, log = 'x')
lines(q_pred)


