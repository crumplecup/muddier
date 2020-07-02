library(data.table)
library(magrittr)

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


# take weighted mean



# order by weighted mean
index <- char_pmfs %>% rownames %>% as.numeric %>% sort
da_mn <- apply(da_df, 2, function (x) weighted.mean(index, x))
dadf <- da_df[, order(da_mn)]
damn <- sort(da_mn)
dmn <- apply(dadf, 2, function (x) weighted.mean(index, x))
# dmn == damn

# convolve sequential flows
it <- array(0, c(nrow(dadf), ncol(dadf)-1))
for (i in 1:(ncol(dadf)-1)) {
  it[ , i] <- convo(dadf[ , i+1], dadf[ , i], index)
}

bit <- convo(dadf[,5], dadf[,4], index)
weighted.mean(index, bit)

itmn <- apply(it, 2, function (x) weighted.mean(index, x))

# fit to exponential distribution
it_exp <- fitdist(itmn, 'exp')
plot(it_exp)
coef(it_exp)

dmn <- damn[-1]
dmn <- dmn[order(itmn)]

plot(log10(damn[-c(1:4)]), itmn[-c(1:3)])
plot(damn[-c(1:4)], itmn[-c(1:3)], log = 'xy', pch = 20, col = get_palette('ocean'),
     xlab = 'deposit age', ylab = 'interarrival time')

# interarrival times vs. deposit ages appear to follow a power law
# representing an evacuation rate of deposits
# the expected interarrival rate is the rate before evacuation

# fit interarrival times to deposit ages to detrend

df <- data.frame(it = log(itmn[-c(1:3)]), age = log(damn[-c(1:4)]))
lit <- lm('it ~ age', data = df)
summary(lit)
plot(lit)

pit <- predict(lit, newdata = data.frame(age = log(damn[-c(1:4)])))

plot(damn[-c(1:4)], itmn[-c(1:3)], log = 'xy', pch = 20, col = get_palette('ocean'),
     xlab = 'deposit age', ylab = 'interarrival time')
lines(damn[-c(1:4)], exp(pit), pch = 20, col = get_palette('forest'), lwd = 2)

# subtract min time from all interarrival times
# fit to deposit ages
# subtract fit times from all interarrival times

df <- data.frame(it = log(itmn[-c(1:3)]),
                 age = log(damn[-c(1:4)]))
lit <- lm('it ~ age', data = df)
summary(lit)

pt <- exp(coef(lit)[1] + log(damn[-c(1:4)]) * coef(lit)[2])
pd <- pt - min(pt)
ip <- itmn[-c(1:3)] - pd
ip[ip <= 0] <- 1
dmn <- damn[-c(1:4)]
plot(dmn, ip, log = 'xy', pch = 20, col = get_palette('ocean'),
     xlab = 'deposit age', ylab = 'interarrival time')
points(dmn[ip == 1], ip[ip == 1], pch = 20, col = get_palette('slate'))
points(dmn[ip == 1], ip[ip == 1], pch = 20, col = get_palette('coral'))

plot(damn[-c(1:4)], itmn[-c(1:3)], log = 'xy', pch = 20, col = get_palette('ocean'),
     xlab = 'deposit age', ylab = 'interarrival time')
lines(damn[-c(1:4)], exp(pit), pch = 20, col = get_palette('forest'), lwd = 3)
lines(damn[-c(1:4)], pd, pch = 20, col = get_palette('coral'), lwd = 3)


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

# ratio of contributing area represented by nodes with debris flows
crk_contr <- br_contr + cd_contr + hf_contr + kn_contr + gr_contr
crk_N <- br_N + cd_N + hf_N + kn_N + gr_N
df_contr <- crk_contr * ncol(dadf) / (crk_N)

ic <- ip * df_contr
ic[ic == df_contr] <- 1

dm500 <- dmn[dmn <= 500]
ic500 <- ic[dmn <= 500]
df500 <- data.frame(it = ic500, age = dm500)
m500 <- lm(it ~ age, data = df500)
summary(m500)

plot(dmn, ic, log = 'xy', pch = 20, col = get_palette('ocean'),
     xlab = 'deposit age', ylab = 'interarrival time per km2')
points(dmn[ic == 1], ip[ip == 1], pch = 20, col = get_palette('slate'))
points(dmn[ic == 1], ip[ip == 1], pch = 20, col = get_palette('coral'))
abline(coef = coef(m500), lwd = 3)
abline(h = mean(ic500), lwd = 2, col = get_palette('charcoal', .8))
text(9500, 125, paste0('mean = ', round(mean(ic500), 2)), cex = .8)

# save objects

setwd('/home/crumplecup/work/muddier/')
usethis::use_data(dadf, overwrite = T)
usethis::use_data(it, overwrite = T)
usethis::use_data(dmn, overwrite = T)
usethis::use_data(ic, overwrite = T)
usethis::use_data(ip, overwrite = T)
usethis::use_data(m500, overwrite = T)


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
da_df_lci <- apply(df, 1, function (x) convo(x+50, rev(df_pmfs_cis[1,]), index))
da_df_uci <- apply(df, 1, function (x) convo(x+50, rev(df_pmfs_cis[3,]), index))
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






