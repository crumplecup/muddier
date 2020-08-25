library(data.table)
library(magrittr)


# sort charcoal pmfs by facies and remove oversamples

index <- char_pmfs %>% rownames %>% as.numeric %>% sort
pmfs <- t(char_pmfs)
ct <- charcoal[ , .N, by = c('facies', 'family')]

df_ct <- ct[ct$facies == 'DF' & ct$N > 2, ]
min_ids <- vector(length = nrow(df_ct), mode = 'character')

for (i in 1:nrow(df_ct)) {
  sub <- charcoal[charcoal$family == unlist(df_ct$family[i]) &
                    charcoal$facies == 'DF', ]
  min_ids[i] <- sub$site_id[sub$mn == min(sub$mn)]
}

rem_ids <- sub$site_id[!(sub$site_id %in% min_ids)]
df <- pmfs[!(charcoal$site_id %in% rem_ids) &
           charcoal$facies == 'DF', ]

ff_ct <- ct[ct$facies == 'FF' & ct$N > 2, ]
min_ids <- vector(length = nrow(ff_ct), mode = 'character')

for (i in 1:nrow(ff_ct)) {
  sub <- charcoal[charcoal$family == unlist(ff_ct$family[i]) &
                    charcoal$facies == 'FF', ]
  min_ids[i] <- sub$site_id[sub$mn == min(sub$mn)]
}

rem_ids <- sub$site_id[!(sub$site_id %in% min_ids)]
ff <- pmfs[!(charcoal$site_id %in% rem_ids) &
             charcoal$facies == 'FF', ]

fg_ct <- ct[ct$facies == 'FG' & ct$N > 2, ]
min_ids <- vector(length = nrow(fg_ct), mode = 'character')

for (i in 1:nrow(fg_ct)) {
  sub <- charcoal[charcoal$family == unlist(fg_ct$family[i]) &
                    charcoal$facies == 'FG', ]
  min_ids[i] <- sub$site_id[sub$mn == min(sub$mn)]
}

rem_ids <- sub$site_id[!(sub$site_id %in% min_ids)]
fg <- pmfs[!(charcoal$site_id %in% rem_ids) &
             charcoal$facies == 'FG', ]


nrow(df)
nrow(ff)
nrow(fg)


# estimate deposit age by subtracting inherited age fit via convolution

index <- char_pmfs %>% rownames %>% as.numeric

begin <- Sys.time()
da_df <- apply(df, 1, function (x) convo(x, rev(df_ph_pmf), index))
end <- Sys.time()
end - begin
nrow(df) / as.numeric(end - begin)

begin <- Sys.time()
da_ff <- apply(ff, 1, function (x) convo(x, rev(ff_ph_pmf), index))
end <- Sys.time()
end - begin
nrow(ff) / as.numeric(end - begin)

begin <- Sys.time()
da_fg <- apply(fg, 1, function (x) convo(x, rev(fg_ph_pmf), index))
end <- Sys.time()
end - begin
nrow(fg) / as.numeric(end - begin)

plot(sort(index), da_ff[,70])

# sum deposit ages and normalize pmf

da_df_pmf <- apply(da_df, 1, function (x) sum(x) / ncol(da_df))
da_ff_pmf <- apply(da_ff, 1, function (x) sum(x) / ncol(da_ff))
da_fg_pmf <- apply(da_fg, 1, function (x) sum(x) / ncol(da_fg))

da_df_cdf <- cumsum(da_df_pmf)
da_ff_cdf <- cumsum(da_ff_pmf)
da_fg_cdf <- cumsum(da_fg_pmf)


# fit cdfs to different distributions and compare GOF

library(fitdistrplus)
descdist(df_cdf %>% as.vector)
df_exp <- fitdist(df_cdf %>% as.vector, 'exp')
df_gam <- fitdist(df_cdf %>% as.vector, 'gamma')
df_bet <- fitdist(df_cdf %>% as.vector, 'beta')
df_wbl <- fitdist(df_cdf %>% as.vector, 'weibull')
df_lnorm <- fitdist(df_cdf %>% as.vector, 'lnorm')


cdfcomp(list(df_exp, df_gam, df_bet, df_wbl, df_lnorm),
        legendtext=c('exp', 'gamma', 'beta', 'weibull', 'lognormal'))
denscomp(list(df_exp, df_gam, df_bet, df_wbl, df_lnorm),
         legendtext=c('exp', 'gamma', 'beta', 'weibull', 'lognormal'))
qqcomp(list(df_exp, df_gam, df_bet, df_wbl, df_lnorm),
       legendtext=c('exp', 'gamma', 'beta', 'weibull', 'lognormal'))
ppcomp(list(df_exp, df_gam, df_bet, df_wbl, df_lnorm),
       legendtext=c('exp', 'gamma', 'beta', 'weibull', 'lognormal'))
gofstat(list(df_exp, df_gam, df_bet, df_wbl, df_lnorm),
        fitnames=c('exp', 'gamma', 'beta', 'weibull', 'lognormal'))

df_exp <- fitdist(df_pmf %>% as.vector, 'exp')

plot(df_exp)
gofstat(df_exp)

# interarrival times for debris flows

# order by weighted mean
index <- char_pmfs %>% rownames %>% as.numeric %>% sort
da_mn <- apply(da_df, 2, function (x) weighted.mean(index, x))
dadf <- da_df[, order(da_mn)]
damn <- sort(da_mn)
dmn <- apply(dadf, 2, function (x) weighted.mean(index, x))
# dmn == damn

# fines
da_mnff <- apply(da_ff, 2, function (x) weighted.mean(index, x))
daff <- da_ff[, order(da_mnff)]
damnff <- sort(da_mnff)
dmnff <- apply(daff, 2, function (x) weighted.mean(index, x))

# gravels
da_mnfg <- apply(da_fg, 2, function (x) weighted.mean(index, x))
dafg <- da_fg[, order(da_mnfg)]
damnfg <- sort(da_mnfg)
dmnfg <- apply(dafg, 2, function (x) weighted.mean(index, x))


# convolve sequential flows
it <- array(0, c(nrow(dadf), ncol(dadf)-1))
for (i in 1:(ncol(dadf)-1)) {
  it[ , i] <- convo(dadf[ , i+1], dadf[ , i], index)
}

itff <- array(0, c(nrow(daff), ncol(daff)-1))
for (i in 1:(ncol(daff)-1)) {
  itff[ , i] <- convo(daff[ , i+1], daff[ , i], index)
}

itfg <- array(0, c(nrow(dafg), ncol(dafg)-1))
for (i in 1:(ncol(dafg)-1)) {
  itfg[ , i] <- convo(dafg[ , i+1], dafg[ , i], index)
}



# take weighted mean
itmn <- apply(it, 2, function (x) weighted.mean(index, x))
itmnff <- apply(itff, 2, function (x) weighted.mean(index, x))
itmnfg <- apply(itfg, 2, function (x) weighted.mean(index, x))

# fit to exponential distribution
it_exp <- fitdist(itmn, 'exp')
plot(it_exp)
coef(it_exp)

dmn <- damn[-1]
dmn <- dmn[order(itmn)]

plot(log10(damn[-c(1:4)]), itmn[-c(1:3)])
plot(damn[-c(1:4)], itmn[-c(1:3)], log = 'xy', pch = 20, col = get_palette('ocean'),
     xlab = 'deposit age', ylab = 'interarrival time')

plot(damnff[-c(1)], itmnff, log = 'xy', pch = 20, col = get_palette('ocean'),
     xlab = 'deposit age', ylab = 'interarrival time')

plot(damnfg[-c(1:2)], itmnfg[-1], log = 'xy', pch = 20, col = get_palette('ocean'),
     xlab = 'deposit age', ylab = 'interarrival time')
points(damn[-c(1:4)], itmn[-c(1:3)], pch = 20, col = get_palette('forest'))
points(damnff[-c(1)], itmnff, pch = 20, col = get_palette('coral'))
legend('topleft', legend = c('debris flows', 'fluvial fines', 'fluvial gravels'),
       fill = get_palette(c('forest', 'coral', 'ocean'), .7))

# interarrival times vs. deposit ages appear to follow a power law
# representing an evacuation rate of deposits
# the expected interarrival rate is the rate before evacuation

# fit interarrival times to deposit ages to detrend

its <- itmn[-c(1:3)]
ages <- damn[-c(1:4)]
its <- its[ages <= 1000]
ages <- ages[ages <= 1000]
df <- data.frame(it = log(its), age = log(ages))

its <- itmnff
ages <- damnff[-c(1)]
its <- its[ages <= 1000]
ages <- ages[ages <= 1000]
dff <- data.frame(it = log(its), age = log(ages))

its <- itmnfg[-c(1)]
ages <- damnfg[-c(1:2)]
its <- its[ages <= 1000]
ages <- ages[ages <= 1000]
dfg <- data.frame(it = log(its), age = log(ages))

lit <- lm('it ~ age', data = df)
flit <- lm('it ~ age', data = dff)
glit <- lm('it ~ age', data = dfg)
summary(lit)
summary(flit)
summary(glit)


pit <- predict(lit, newdata = data.frame(age = log(damn[-c(1:4)])))
pitf <- predict(flit, newdata = data.frame(age = log(damn[-c(1)])))
pitg <- predict(glit, newdata = data.frame(age = log(damn[-c(1:2)])))



png('it1.png', height = 15, width = 21, units = 'cm', res = 300)
plot(damnfg[-c(1:2)], itmnfg[-1], log = 'xy', pch = 20, col = get_palette('ocean'),
     xlab = 'deposit age', ylab = 'interarrival time')
points(damn[-c(1:4)], itmn[-c(1:3)], pch = 20, col = get_palette('forest'))
points(damnff[-c(1)], itmnff, pch = 20, col = get_palette('coral'))
legend('topleft', legend = c('debris flows', 'fluvial fines', 'fluvial gravels'),
       fill = get_palette(c('forest', 'coral', 'ocean'), .7))
lines(damn[-c(1:4)], exp(pit), pch = 20, col = get_palette('forest'), lwd = 2, lty = 3)
lines(damn[-c(1)], exp(pitf), pch = 20, col = get_palette('coral'), lwd = 2, lty = 3)
lines(damn[-c(1:2)], exp(pitg), pch = 20, col = get_palette('ocean'), lwd = 2, lty = 3)
dev.off()

# subtract min time from all interarrival times
# fit to deposit ages
# subtract fit times from all interarrival times


pt <- exp(coef(lit)[1] + log(damn[-c(1:4)]) * coef(lit)[2])
pd <- pt - min(pt)
ip <- itmn[-c(1:3)] - pd
ip[ip <= 0] <- 1
dmn <- damn[-c(1:4)]

fpt <- exp(coef(flit)[1] + log(damnff[-c(1)]) * coef(flit)[2])
fpd <- fpt - min(fpt)
fip <- itmnff - fpd
fip[fip <= 0] <- 1
fdmn <- damnff[-c(1)]

gpt <- exp(coef(glit)[1] + log(damnfg[-c(1:2)]) * coef(glit)[2])
gpd <- gpt - min(gpt)
gip <- itmnfg[-c(1)] - gpd
gip[gip <= 0] <- 1
gdmn <- damnfg[-c(1:2)]

plot(dmn, ip, log = 'xy', pch = 20, col = get_palette('ocean'),
     xlab = 'deposit age', ylab = 'interarrival time')
points(dmn[ip == 1], ip[ip == 1], pch = 20, col = get_palette('slate'))
points(dmn[ip == 1], ip[ip == 1], pch = 20, col = get_palette('coral'))

plot(damn[-c(1:4)], itmn[-c(1:3)], log = 'xy', pch = 20, col = get_palette('ocean'),
     xlab = 'deposit age', ylab = 'interarrival time')
lines(damn[-c(1:4)], exp(pit), pch = 20, col = get_palette('forest'), lwd = 3)
lines(damn[-c(1:4)], pd, pch = 20, col = get_palette('coral'), lwd = 3)


# temporal scoping by floating background window

# i try four different types of background window for debris flows
# finding i like the 'log' method best, i fit fluvial fines and gravels with it

# organize observations into windows of time
bw <- back_window(itmn[-c(1:3)], dmn)
bwl <- back_window(itmn[-c(1:3)], dmn, type = 'log')
bwfl <- back_window(itmn[-c(1:3)], dmn, type = 'float_log')
bwd <- back_window(itmn[-c(1:3)], dmn, type = 'discrete')

# fit interarrival time (it) to deposit age and predict expected it for each window type
bwn <- lapply(bwfl, nrow) %>% unlist
bwf <- lapply(bw,
              function(x) lm('it ~ age',
                             data = data.frame(it = log(x[, 1]), age = log(x[, 2]))))
bwp <- mapply(function(x, y) exp(predict(x, newdata = data.frame(age = log(y)))), bwf, dmn)

bwl <- list(bwl[[2]], bwl[[3]], bwl[[4]])
bwlf <- lapply(bwl,
              function(x) lm('it ~ age',
                             data = data.frame(it = log(x[, 1]), age = log(x[, 2]))))
bwlp <- mapply(function(x, y) predict(x, newdata = data.frame(age = log(y[,2]))), bwlf, bwl) %>% unlist
bwlp <- c(bwlp, predict(bwlf[[3]], newdata = data.frame(age = log(dmn[length(dmn)]))))

bwflf <- lapply(bwfl,
              function(x) lm('it ~ age',
                             data = data.frame(it = log(x[, 1]), age = log(x[, 2]))))
bwflp <- mapply(function(x, y) exp(predict(x, newdata = data.frame(age = log(y)))), bwflf, dmn)


bwd <- list(bwd[[1]], bwd[[2]])
its <- itmn[-c(1:3)]
bwd[[3]] <- matrix(c(its[dmn >= 1000], as.numeric(dmn[dmn >= 1000])), ncol = 2)
bwdf <- lapply(bwd,
              function(x) lm('it ~ age',
                             data = data.frame(it = log(x[, 1]), age = log(x[, 2]))))
bwdp <- mapply(function(x, y) predict(x, newdata = data.frame(age = log(y[,2]))), bwdf, bwd) %>% unlist

# fluvial fines, log method
ffbwl <- back_window(itmnff, dmnff[-1], type = 'log')
ffbwl <- list(ffbwl[[2]], ffbwl[[3]], ffbwl[[4]])
ffbwf <- lapply(ffbwl,
              function(x) lm('it ~ age',
                             data = data.frame(it = log(x[, 1]), age = log(x[, 2]))))
ffbwp <- mapply(function(x, y) predict(x, newdata = data.frame(age = log(y[,2]))), ffbwf, ffbwl) %>% unlist
ffbwp <- c(ffbwp, predict(ffbwf[[3]], newdata = data.frame(age = log(dmnff[(length(dmnff)-1):length(dmnff)]))))

plot(dmnff[-1], itmnff, log = 'xy')
points(dmnff[-1], exp(ffbwp), pch = 20, col = 'slateblue')

# fluvial gravels, log method
fgbwl <- back_window(itmnfg, dmnfg[-1], type = 'log')
fgbwl <- list(fgbwl[[2]], fgbwl[[3]], fgbwl[[4]])
fgbwf <- lapply(fgbwl,
                function(x) lm('it ~ age',
                               data = data.frame(it = log(x[, 1]), age = log(x[, 2]))))
fgbwp <- mapply(function(x, y) predict(x, newdata = data.frame(age = log(y[,2]))), fgbwf, fgbwl) %>% unlist
fgbwp <- c(
#  predict(fgbwf[[1]], newdata = data.frame(age = log(dmnfg[1]))),
  fgbwp,
  predict(fgbwf[[3]], newdata = data.frame(age = log(dmnfg[(length(dmnfg)-3):length(dmnfg)]))))
plot(dmnfg[-c(1:2)], itmnfg[-1], log = 'xy')
points(dmnfg[-c(1:2)], exp(fgbwp), pch = 20, col = 'slateblue')

png('da_it.png', height = 17, width = 21, units = 'cm', res = 300)
plot(dmnfg[-c(1:2)], itmnfg[-1], log = 'xy', pch = 20, col = get_palette('ocean'),
     xlab = 'deposit age', ylab = 'interarrival time')
lines(dmnfg[-c(1:2)], exp(fgbwp), lwd = 2, lty = 3, col = get_palette('ocean'))
points(dmnff[-1], itmnff, pch = 20, col = get_palette('coral'))
lines(dmnff[-1], exp(ffbwp), lwd = 2, lty = 3, col = get_palette('coral'))
points(dmn, itmn[-c(1:3)], pch = 20, col = get_palette('forest'))
lines(dmn, exp(bwlp), lwd = 2, lty = 3, col = get_palette('forest'))
legend('topleft', legend = c('debris flows', 'fluvial fines', 'fluvial gravels'),
       fill = get_palette(c('forest', 'coral', 'ocean'), .7))
dev.off()


ldf_pd <- exp(bwlp) - exp(min(bwlp))
ldf_ip <- itmn[-c(1:3)] - ldf_pd
mean(ldf_ip[dmn <= 100])
mean(ldf_ip[dmn > 100 & dmn <= 1000])
mean(ldf_ip[dmn > 1000])
rt_df <- mean(ldf_ip[dmn <= 1000])
ex_df <- dmn / rt_df


ct_df <- 1:length(damn[-1])
ct_ff <- 1:length(dmnff[-1])
ct_fg <- 1:length(dmnfg[-1])
plot(dmnfg[-c(1:2)], dmnfg[-c(1:2)] / ct_fg[-1] * samp_contr, pch = 20,
     col = get_palette('ocean'), log = 'xy',
     xlab = 'deposit age', ylab = 'mean interarrival time per km2 (age <= x)')
points(damn[-1], damn[-1] / ct_df * samp_contr, pch = 20, col = get_palette('forest'))
points(dmnff[-1], dmnff[-1] / ct_ff * samp_contr, pch = 20, col = get_palette('coral'))

png('da_mnit.png', height = 17, width = 21, units = 'cm', res = 300)
plot(dmnfg[-c(1:2)], dmnfg[-c(1:2)] / ct_fg[-1], pch = 20,
     col = get_palette('ocean'), log = 'xy',
     xlab = 'deposit age', ylab = 'mean interarrival time (age <= x)')
points(damn[-1], damn[-1] / ct_df, pch = 20, col = get_palette('forest'))
points(dmnff[-1], dmnff[-1] / ct_ff, pch = 20, col = get_palette('coral'))
legend('topleft', legend = c('debris flows', 'fluvial fines', 'fluvial gravels'),
       fill = get_palette(c('forest', 'coral', 'ocean'), .7))
dev.off()


# the rate at which the interarrival period increases represents
# the removal rate of the stream
# we can represent the removal rate using the hazard functin of deposit ages
# the hazard function is the pmf divided by the exceedance function

dadf_pmf <- apply(dadf, 1, function(x) sum(x) / length(x))
dadf_exc <- to_exceed(dadf_pmf)
dadf_hzd <- dadf_pmf / dadf_exc
# dadf_hzd[dadf_hzd == max(dadf_hzd)] <- 3
# dadf_hzd[dadf_hzd > .1] <- .1

daff_pmf <- apply(daff, 1, function(x) sum(x) / length(x))
daff_exc <- to_exceed(daff_pmf)
daff_hzd <- daff_pmf / daff_exc
# daff_hzd[daff_hzd < 0] <- 0

dafg_pmf <- apply(dafg, 1, function(x) sum(x) / length(x))
dafg_exc <- to_exceed(dafg_pmf)
dafg_hzd <- dafg_pmf / dafg_exc


index <- char_pmfs %>% rownames %>% as.numeric %>% sort
options(scipen = 5)
setwd('/home/crumplecup/work/')
png('dadf_hzd.png', height = 17, width = 29, units = 'cm', res = 300)
plot(index, dadf_hzd, log = 'xy', pch = 20, col = get_palette('crimson'),
     xlab = 'Debris Flow Deposit Age', ylab = 'Hazard Function')
dev.off()

png('daff_hzd.png', height = 17, width = 29, units = 'cm', res = 300)
plot(index, daff_hzd, log = 'xy', pch = 20, col = get_palette('crimson'),
     xlab = 'Fluvial Fines Deposit Age', ylab = 'Hazard Function')
dev.off()

png('dafg_hzd.png', height = 17, width = 29, units = 'cm', res = 300)
plot(index, dafg_hzd, log = 'xy', pch = 20, col = get_palette('crimson'),
     xlab = 'Fluvial Gravels Deposit Age', ylab = 'Hazard Function')
dev.off()

pal <- get_palette(c('crimson', 'charcoal', 'ocean'), .8)
png('da_hzd.png', height = 17, width = 29, units = 'cm', res = 300)
plot(index, dafg_hzd, log = 'xy', type = 'l', lwd = 2, col = pal[3],
     xlab = 'Deposit Age', ylab = 'Hazard Function')
lines(index, daff_hzd, lty = 2, lwd = 2, col = pal[2])
lines(index, dadf_hzd, pch = 20, lwd = 2, col = pal[1])
legend('topleft', legend = c('Debris Flows', 'Fluvial Fines', 'Fluvial Gravels'),
       fill = get_palette(c('crimson', 'charcoal', 'ocean'), .8))
dev.off()


# fit an exponential distribution to mean estimated interarrival times
# that rate is output rate, because the whole record would be present if only input
# figure out which 5-year charcoal age increment each deposit age would fall under
# round the mean estimated age to the nearest 5-year increment

dadf_mn <- apply(dadf, 2, function(x) weighted.mean(index, x))
daff_mn <- apply(daff, 2, function(x) weighted.mean(index, x))
dafg_mn <- apply(dafg, 2, function(x) weighted.mean(index, x))


# events stores the number of events in each 5-year increment
# expand that data into a vector of ages, one age for each event
# to estimate the time between events

index <- char_pmfs %>% rownames %>% as.numeric %>% sort
dadf_it <- arrivals(dadf_mn, index)
daff_it <- arrivals(daff_mn, index)
dafg_it <- arrivals(dafg_mn, index)

# pmf, cdf of interarrival time

dadf_it_cdf <- emp_cdf(dadf_it)
daff_it_cdf <- emp_cdf(daff_it)
dafg_it_cdf <- emp_cdf(dafg_it)

dadf_it_mod <- lm('it ~ age', data = data.frame(it = dadf_it_cdf[,2], age = log(dadf_it_cdf[,1] + .0001)))
daff_it_mod <- lm('it ~ age', data = data.frame(it = daff_it_cdf[,2], age = log(daff_it_cdf[,1] + .0001)))
dafg_it_mod <- lm('it ~ age', data = data.frame(it = dafg_it_cdf[,2], age = log(dafg_it_cdf[,1] + .0001)))

dadf_it_prd <- predict(dadf_it_mod, newdata = data.frame(age = log(dadf_it_cdf[,1]+.0001)))
daff_it_prd <- predict(daff_it_mod, newdata = data.frame(age = log(daff_it_cdf[,1]+.0001)))
dafg_it_prd <- predict(dafg_it_mod, newdata = data.frame(age = log(dafg_it_cdf[,1]+.0001)))
summary(dadf_it_mod)
summary(daff_it_mod)
summary(dafg_it_mod)

plot(dadf_it_cdf[,1], dadf_it_prd / max(dadf_it_prd),
     type = 'l', ylim = c(0,1), lwd = 2, col = pal[1], log = 'x',
     xlab = 'Interarrival Time', ylab = 'CDF')
points(dadf_it_cdf[,1], dadf_it_cdf[,2]+.001, pch = 20, col = pal[1])
lines(daff_it_cdf[,1], daff_it_prd / max(daff_it_prd), lwd = 2, col = pal[2])
points(daff_it_cdf[,1], daff_it_cdf[,2]+.001, pch = 20, col = pal[2])
lines(dafg_it_cdf[,1], dafg_it_prd / max(dafg_it_prd), lwd = 2, col = pal[3])
points(dafg_it_cdf[,1], dafg_it_cdf[,2]+.001, pch = 20, col = pal[3])

dadf_it_pmf <- to_pmf(dadf_it_cdf[,2])
daff_it_pmf <- to_pmf(daff_it_cdf[,2])
dafg_it_pmf <- to_pmf(dafg_it_cdf[,2])

dadf_it_exp <- fitdist(as.numeric(dadf_it), 'exp')
dadf_it_cdf_fit <- 1 - exp(-dadf_it_exp$estimate * dadf_it_cdf[,1])
plot(dadf_it_cdf[,1], dadf_it_cdf_fit, log = 'x')

# fit an exponential distribution to arrivals < 1000
# older arrivals are disappearing from the deposit record, deflating the rate

df_ages <- dadf_mn[-1]
dadf_it_exp_yng <- fitdist(as.numeric(dadf_it[df_ages <= 1000]), 'exp')
# results in mean expected interarrival time of 8.18 years
# shows that older deposits are disappearing from the record
# the rate at which deposits are disappearing
# describes the size of the initial pool

plot(df_ages, dadf_it, log = 'xy')



# fit exp dist to deposit ages
# find the proportion of deposits remaining using cdf of ages
# weight the interarrival times by the proportion of remaining deposits

dadf_age <- ages(dadf_mn, index) %>% as.numeric
dadf_exp <- fitdist(dadf_age/10, 'exp')
dadf_age_prd <- 1 - exp(-(dadf_exp$estimate/10) * (dadf_age+.15))
dadf_itc <- dadf_it * (1 - dadf_age_prd[-1])
dadf_itc_exp <- fitdist(dadf_itc %>% as.numeric, 'exp')

daff_age <- ages(daff_mn, index) %>% as.numeric
daff_exp <- fitdist(daff_age/10, 'exp')
daff_age_prd <- 1 - exp(-(daff_exp$estimate/10) * (daff_age+.15))
daff_itc <- daff_it * (1 - daff_age_prd[-1])
daff_itc_exp <- fitdist(daff_itc %>% as.numeric, 'exp')

dafg_age <- ages(dafg_mn, index) %>% as.numeric
dafg_exp <- fitdist(dafg_age/10, 'exp')
dafg_age_prd <- 1 - exp(-(dafg_exp$estimate/10) * (dafg_age+.15))
dafg_itc <- dafg_it * (1 - dafg_age_prd[-1])
dafg_itc_exp <- fitdist(dafg_itc %>% as.numeric, 'exp')

plot(dadf_age[-1], dadf_it, log = 'xy',
     xlab = 'Deposit Age', ylab = 'Interarrival Time')

plot(dadf_age[-1], dadf_itc, log = 'x')
plot(daff_age[-1], daff_itc, log = 'x')
plot(dafg_age[-1], dafg_itc, log = 'x')

png('age_itc_ln.png', height = 17, width = 33, units = 'cm', res = 300)
pal <- get_palette(c('crimson', 'charcoal', 'ocean'), .4)
plot(daff_age[-1], daff_itc+1, pch = 20, col = pal[2], log = 'xy',
     xlab = 'Deposit Age', ylab = 'Interarrival Times + 1 (discounted)')
points(dadf_age[-1], dadf_itc+1, pch = 20, col = pal[1])
points(dafg_age[-1], dafg_itc+1, pch = 20, col = pal[3])
dev.off()

png('age_itc.png', height = 17, width = 33, units = 'cm', res = 300)
pal <- get_palette(c('crimson', 'charcoal', 'ocean'), .4)
plot(daff_age[-1], daff_itc, pch = 20, col = pal[2], log = 'x',
     xlab = 'Deposit Age', ylab = 'Interarrival Times (discounted)')
points(dadf_age[-1], dadf_itc, pch = 20, col = pal[1])
points(dafg_age[-1], dafg_itc, pch = 20, col = pal[3])
dev.off()


png('age_it.png', height = 17, width = 33, units = 'cm', res = 300)
pal <- get_palette(c('crimson', 'charcoal', 'ocean'), .4)
plot(dadf_age[-1], dadf_it+1, log = 'xy', pch = 20, col = pal[1],
     xlab = 'Deposit Age', ylab = 'Interarrival Time (+1)')
points(daff_age[-1], daff_it+1, pch = 20, col = pal[2])
points(dafg_age[-1], dafg_it+1, pch = 20, col = pal[3])
legend('topleft', legend = c('Debris Flows', 'Fluvial Fines', 'Fluvial Gravels'),
       fill = get_palette(c('crimson', 'charcoal', 'ocean'), .8))
dev.off()

plot(index, dadf_exp_prd, log = 'x')
points(index, cumsum(dadf_pmf), pch = 20, col = pal[1])
lines(index, dadf_cdf, lwd = 2, col = pal[2])
lines(index, dadf_cdf_prd, lwd = 3, col = pal[3])
points(index, dadf_mns)

shp <- dadf_wbl$estimate[1]
scl <- dadf_wbl$estimate[2]

dadf_wbl_prd <- (shp / scl) *
  ((index + .001) / scl)^(shp - 1) *
  exp(-((index + .001) / scl)^shp)



dadf_it_pmf_mod <- lm('it ~ age', data = data.frame(it = log(dadf_it_pmf+.001), age = log(dadf_it_cdf[,1]+.001)))
summary(dadf_it_pmf_mod)
dadf_it_pmf_prd <- predict(dadf_it_pmf_mod, newdata = data.frame(age = log(dadf_it_cdf[,1]+.001)))
daff_it_pmf_prd <- to_pmf(daff_it_prd / max(daff_it_prd))
dafg_it_pmf_prd <- to_pmf(dafg_it_prd / max(dafg_it_prd))

plot(dadf_it_cdf[,1], dadf_it_pmf, log = 'x', pch = 20, col = pal[1])
lines(dadf_it_cdf[,1], exp(dadf_it_pmf_prd), lwd = 2, col = pal[1])

plot(dadf_it_cdf[,1]+1, dadf_it_prd / max(dadf_it_prd), pch = 20, col = 'slateblue', log = 'x')
points(dadf_it_cdf[,1]+1, dadf_it_cdf[,2]+.001, pch = 20, col = get_palette('crimson'))

dadf_exp <- lm('it ~ age', data = data.frame(it = log(dadf_it), age = log(dadf_mn[-1])))
summary(dadf_exp)
dadf_prd <- predict(dadf_exp, newdata = data.frame(age = log(index)))

plot(dadf_mn[-1], dadf_it+1, log = 'xy')
plot(dadf_mn[-1], log = 'y')

# scope of interarrival period

# change in contributing area over each study area

br_contr <- max(creeks$contr_area[creeks$creek_name == 'bear'])
cd_contr <- max(creeks$contr_area[creeks$creek_name == 'cedar'])
hf_contr <- max(creeks$contr_area[creeks$creek_name == 'hoffman'])
kn_contr <- max(creeks$contr_area[creeks$creek_name == 'knowles'])
gr_contr <- 4.732 ## from GRC_SedVol_new_20170315

# number of nodes per study area
br_N <- nrow(creeks[creeks$creek_name == 'bear', ])
cd_N <- nrow(creeks[creeks$creek_name == 'cedar', ])
hf_N <- nrow(creeks[creeks$creek_name == 'hoffman', ])
kn_N <- nrow(creeks[creeks$creek_name == 'knowles', ])
gr_N <- round(kn_N * gr_contr / kn_contr)

# sample sites per study area

br_fcn <- charcoal$facies[grep('BC', charcoal$family)]
br_dfn <-length(grep('DF', br_fcn))
br_ffn <- length(grep('FF', br_fcn))
br_fgn <- length(grep('FG', br_fcn))

cd_fcn <- charcoal$facies[grep('CC', charcoal$family)]
cd_dfn <-length(grep('DF', cd_fcn))
cd_ffn <- length(grep('FF', cd_fcn))
cd_fgn <- length(grep('FG', cd_fcn))

gr_fcn <- c(
  charcoal$facies[grep('GR', charcoal$family)],
  charcoal$facies[grep('T8', charcoal$family)])
gr_dfn <-length(grep('DF', gr_fcn))
gr_ffn <- length(grep('FF', gr_fcn))
gr_fgn <- length(grep('FG', gr_fcn))

kn_fcn <- c(charcoal$facies[grep('UK', charcoal$family)],
            charcoal$facies[grep('LK', charcoal$family)],
            charcoal$facies[grep('Dam', charcoal$family)],
            charcoal$facies[grep('DF', charcoal$family)])
kn_dfn <-length(grep('DF', kn_fcn))
kn_ffn <- length(grep('FF', kn_fcn))
kn_fgn <- length(grep('FG', kn_fcn))

charcoal$site_id[!charcoal$site_id %in% charcoal$site_id[c(
  grep('UK', charcoal$family),
  grep('LK', charcoal$family),
  grep('Dam', charcoal$family),
  grep('DFK', charcoal$family),
  grep('GR', charcoal$family),
  grep('T8', charcoal$family),
  grep('BC', charcoal$family),
  grep('CC', charcoal$family)
)]]

# ratio of contributing area represented by nodes with debris flows
crk_contr <- br_contr + cd_contr + hf_contr + kn_contr + gr_contr
crk_N <- br_N + cd_N + hf_N + kn_N + gr_N
df_contr <- crk_contr * ncol(dadf) / (crk_N)
ff_contr <- crk_contr * ncol(daff) / (crk_N)
fg_contr <- crk_contr * ncol(dafg) / (crk_N)
samp_contr <- crk_contr * (ncol(dadf) + ncol(daff) + ncol(dafg))/ crk_N

ic <- ip * df_contr
ic[ic == df_contr] <- 1
fic <- fip * ff_contr
fic[fic == ff_contr] <- 1.2
gic <- gip * fg_contr
gic[gic == fg_contr] <- 1.5


png('cor_it.png', height = 15, width = 21, units = 'cm', res = 300)
plot(dmn, ic, log = 'xy', pch = 20, col = get_palette('forest'),
     xlab = 'deposit age', ylab = 'interarrival time per km2')
abline(h = mean(ic[dmn <= 1000]), lwd = 2, col = get_palette('forest', .8), lty = 3)
text(10700, 110, paste0('mean = ', round(mean(ic[dmn <= 1000]), 2)), cex = .8,
     col = get_palette('forest', 1))
points(fdmn, fic, pch = 20, col = get_palette('coral'))
abline(h = mean(fic[fdmn <= 1000]), lwd = 2, col = get_palette('coral', .8), lty = 3)
text(10700, 16, paste0('mean = ', round(mean(fic[fdmn <= 1000]), 2)), cex = .8,
     col = get_palette('coral', 1))
points(gdmn, gic, pch = 20, col = get_palette('ocean'))
abline(h = mean(gic[gdmn <= 1000]), lwd = 2, col = get_palette('ocean', .8), lty = 3)
text(10700, 45, paste0('mean = ', round(mean(gic[gdmn <= 1000]), 2)), cex = .8,
     col = get_palette('ocean', 1))
legend('topleft', legend = c('debris flows', 'fluvial fines', 'fluvial gravels'),
       fill = get_palette(c('forest', 'coral', 'ocean'), .7))
dev.off()

# save objects

setwd('/home/crumplecup/work/muddier/')
usethis::use_data(dadf, overwrite = T)
usethis::use_data(it, overwrite = T)
usethis::use_data(dmn, overwrite = T)
usethis::use_data(ic, overwrite = T)
usethis::use_data(ip, overwrite = T)

usethis::use_data(daff, overwrite = T)
usethis::use_data(itff, overwrite = T)
usethis::use_data(fdmn, overwrite = T)
usethis::use_data(fic, overwrite = T)
usethis::use_data(fip, overwrite = T)

usethis::use_data(dafg, overwrite = T)
usethis::use_data(itfg, overwrite = T)
usethis::use_data(gdmn, overwrite = T)
usethis::use_data(gic, overwrite = T)
usethis::use_data(gip, overwrite = T)


# fire interarrival times

index <- char_pmfs %>% rownames %>% as.numeric
ca <- apply(char_pmfs[,-1], 2, function(x) weighted.mean(index, x))

# convolve sequential flows
mat <- as.matrix(char_pmfs)
fit <- array(0, c(nrow(mat), ncol(mat)-1))
for (i in 1:(ncol(mat)-1)) {
  fit[ , i] <- convo(mat[ , i+1], mat[ , i], index)
}

index <- index %>% sort
fitmn <- apply(fit, 2, function(x) weighted.mean(index, x))

fdf <- data.frame(it = log(fitmn), age = log(ca))
mit <- lm('it ~ age', data = fdf)
summary(mit)

mp <- predict(mit, newdata = data.frame(age = log(ca)))

mpt <- exp(coef(mit)[1] + log(ca) * coef(mit)[2])
mpd <- mpt - min(mpt)
mip <- fitmn - mpd
mip[mip <= 1] <- 1

png('fire_it.png', height = 17, width = 21, units = 'cm', res = 300)
plot(ca, fitmn, log = 'xy', pch = 20, col = get_palette('crimson'),
     xlab = 'charcoal age', ylab = 'interarrival time',
     main = 'fire interarrival times')
dev.off()

plot(ca, mip, log = 'xy', pch = 20, col = get_palette('crimson'))

plot(density(dmn), col = get_palette('crimson', .9))
lines(density(ca))

# older plots

pal <- get_palette(c('ocean', 'sky', 'gold', 'coral', 'forest', 'charcoal'), .7)
plot(sort(index), da_df_cdf, type = 'l', lwd = 2, col = pal[5], ylim = c(0,1))
lines(sort(index), da_ff_cdf, lwd = 2, col = pal[1])
lines(sort(index), da_fg_cdf, lwd = 2, col = pal[3])

da_df_exp <- fitdistrplus::fitdist(da_df_cdf, 'exp')
da_df_wbl <- fitdistrplus::fitdist(da_df_cdf, 'weibull')
lines(sort(index), da_df_exp$data)
plot(da_df_exp)
fitdistrplus::descdist(da_df_cdf)

bit <- da_df_exp$estimate
bot <- -bit * sort(index)
ebot <- exp(bot)

png('dep_ages1.png', height = 5, width = 8, units = 'in', res = 300)
plot(sort(index)[-1], da_df_pmf[-1], type = 'l', lwd = 3, col = pal[6],
     log = 'x', xlab = 'Deposit Age', ylab = 'Probability of Age(x)',
     main = 'Distribution of Non-Zero Deposit Ages')
lines(sort(index)[-1], da_ff_pmf[-1], lwd = 3, col = pal[2])
lines(sort(index)[-1], da_fg_pmf[-1], lwd = 3, col = pal[4])
legend('topright', legend = c('Debris Flows', 'Fluvial Fines', 'Fluvial Gravels'),
       fill = c(pal[c(6,2,4)]))
dev.off()

png('inher_age_correction.png', height = 5, width = 7, units = 'in', res = 300)
plot(index, ff[60, ], xlim = c(0, 6000), type = 'l', lwd =2, col = pal[1],
     ylim = c(0, .05), xlab = 'Deposit Age', ylab = 'Probability of Age(x)')
lines(sort(index), da_ff[, 60], lwd = 2, col = pal[2])
lines(index, fg[85, ], lwd = 2, col = pal[3])
lines(sort(index), da_fg[, 85], lwd = 2, col = pal[4])
lines(index, df[150, ], lwd = 2, col = pal[5])
lines(sort(index), da_df[, 150], lwd = 2, col = pal[6])
legend('topright', legend = c('sample age FF', 'deposit age FF',
                              'sample age FG', 'deposit age FG',
                              'sample age DF', 'deposit age DF'),
       fill = pal)
dev.off()


begin <- Sys.time()
da_df_lci <- apply(df, 1, function (x) convo(x, rev(df_pmfs_cis[1,]), index))
da_df_uci <- apply(df, 1, function (x) convo(x, rev(df_pmfs_cis[3,]), index))
end <- Sys.time()
end - begin
nrow(df) / as.numeric(end - begin)

begin <- Sys.time()
da_ff_lci <- apply(ff, 1, function (x) convo(x, rev(ff_pmfs_cis[, 1]), index))
da_ff_uci <- apply(ff, 1, function (x) convo(x, rev(ff_pmfs_cis[, 3]), index))
end <- Sys.time()
end - begin
nrow(ff) / as.numeric(end - begin)

begin <- Sys.time()
da_fg_lci <- apply(fg, 1, function (x) convo(x, rev(fg_pmfs_cis[1,]), index))
da_fg_uci <- apply(fg, 1, function (x) convo(x, rev(fg_pmfs_cis[3,]), index))
end <- Sys.time()
end - begin
nrow(fg) / as.numeric(end - begin)

dfdap_uci <- da_df_uci
dfdap_uci[dfdap_uci > 1] <- 1
dfdap_uci <- apply(dfdap_uci, 1, function (x) sum(x) / ncol(dfdap_uci))

ffdap_uci <- da_ff_uci
ffdap_uci[ffdap_uci > 1] <- 1
ffdap_uci <- apply(ffdap_uci, 1, function (x) sum(x) / ncol(ffdap_uci))

fgdap_uci <- da_fg_uci
fgdap_uci[fgdap_uci > 1] <- 1
fgdap_uci <- apply(fgdap_uci, 1, function (x) sum(x) / ncol(fgdap_uci))

dfdap_lci <- da_df_lci
dfdap_lci[dfdap_lci > 1] <- 1
dfdap_lci <- apply(dfdap_lci, 1, function (x) sum(x) / ncol(dfdap_lci))

ffdap_lci <- da_ff_lci
ffdap_lci[ffdap_lci > 1] <- 1
ffdap_lci <- apply(ffdap_lci, 1, function (x) sum(x) / ncol(ffdap_lci))

fgdap_lci <- da_fg_lci
fgdap_lci[fgdap_lci > 1] <- 1
fgdap_lci <- apply(fgdap_lci, 1, function (x) sum(x) / ncol(fgdap_lci))


years <- sort(index) + 50
png('dep_ages_df1.png', height = 5, width = 6, units = 'in', res = 300)
plot(years[-1], dfdap_uci[-1], type = 'l', lwd = 2, col = pal[5], log = 'x', # xlim = c(0, 15000),
     xlab = 'Deposit Age', ylab = 'Probability of Age(x)',
     main = 'Uncertainty for Non-Zero Debris-Flow Deposit Ages')
lines(years[-1], da_df_pmf[-1], lwd = 2, col = pal[6])
lines(years[-1], dfdap_lci[-1], lwd = 2, col = pal[5])
legend('topright', legend = c('Debris Flows', '95% CIs'),
       fill = pal[c(6,5)])
dev.off()

png('zero_ages.png', height = 5, width = 6, units = 'in', res = 300)
plot(years[1], da_df_pmf_lci[1], bg = pal[5], pch = 25, ylim = c(0,1),
     xlab = 'Deposit Age', ylab = 'Probability of Age(x)', col = pal[6],
     main = 'Uncertainty for Zero-Age Deposits')
points(years[1], da_df_pmf[1], pch = 21, bg = pal[6])
points(years[1], da_df_pmf_uci[1], pch = 24, bg = pal[5], col = pal[6])

points(years[1] - .07, da_ff_pmf[1], pch = 21, bg = pal[2])
points(years[1] - .07, ffdap_lci[1], pch = 24, bg = pal[5], col = pal[2])
points(years[1] - .07, ffdap_uci[1], pch = 25, bg = pal[5], col = pal[2])

points(years[1] + .07, da_fg_pmf[1], pch = 21, bg = pal[4])
points(years[1] + .07, fgdap_lci[1], pch = 24, bg = pal[5], col = pal[4])
points(years[1] + .07, fgdap_uci[1], pch = 25, bg = pal[5], col = pal[4])
legend('topleft', legend = c(
  'Fluvial Fines', 'Fluvial Gravel', 'Debris Flows', '95% CIs'),
  fill = pal[c(2,4,6,5)])
dev.off()

trc_df_lci <- da_df_cdf_lci
trc_df_lci[trc_df_lci > 1] <- 1

plot(years, trc_df_lci, type = 'l', lwd = 2, col = pal[5],
     xlim = c(0,15000), ylim = c(0, 1))
lines(years, da_df_cdf, lwd = 2, col = pal[6])
lines(years, da_df_cdf_uci, lwd = 2, col = pal[5])


vec <- 1:10000
yv <- .002 * exp(-.002*vec)
plot(vec, yv)
weighted.mean(vec, yv)


# sample ages come from the lab age as years older than 1950
# sample ages cannot be younger than 2000
# we can convert sample ages to year 2020
# but no sample age can occur from 2000 to 2020

# when calculating inherited age
# no sample can be younger than zero age
# which is year 2000, not 1950
# so if sample age +50 < 0, impossible






