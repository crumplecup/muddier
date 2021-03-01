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
  itmn[i-1] <- dfmn[i] - dfmn[i-1]
}

ffmn <- apply(ff, 1, function (x) weighted.mean(rev(index), x))
ffitmn <- 0
for (i in 2:length(ffmn)) {
  ffitmn[i-1] <- ffmn[i] - ffmn[i-1]
}

fgmn <- apply(fg, 1, function (x) weighted.mean(rev(index), x))
fgitmn <- 0
for (i in 2:length(fgmn)) {
  fgitmn[i-1] <- fgmn[i] - fgmn[i-1]
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



library(lubridate)

flow <- fread('/media/crumplecup/catacomb/research/siuslaw_flow.csv')
flow$dt <- mdy_hm(flow$datetime)
flow$yr <- year(flow$dt)
flow$mt <- month(flow$dt)
flow[,max(discharge_cfs), by = yr] %>% plot
flow$dr <- 0
flow$dr[-1] <- flow$dt[-1] - flow$dt[-length(flow$dt)]
flow$dr[1] <- flow$dr[2]
flow$dr[flow$dr > 1800] <- 1800
flow$dq <- 0
flow$dq[-1] <- flow$discharge_cfs[-1] - flow$discharge_cfs[-length(flow$discharge_cfs)]



q <- flow$discharge_cfs %>% sort
plot(emp_cdf(q), pch = 20, col = pal[2],
     xlab = 'discharge in cfs', ylab = 'CDF')


length(tdif[tdif == 900]) / length(tdif) * 900 +
length(tdif[tdif == 1800]) / length(tdif) * 1800

plot(flow$dt[flow$yr == 2020 & flow$mt %in% c(12)],
     flow$dq[flow$yr == 2020 & flow$mt %in% c(12)], type = 'l')


flow$dq[flow$yr == 2020 & flow$mt %in% c(12)]

is_peak <- function(vec) {
  vec[1] > 0 & vec[2] < 0
}

dec20 <- flow[flow$yr == 2020 & flow$mt == 12, ]
plot(dec20$discharge_cfs, type = 'l')
match(max(dec20$discharge_cfs), dec20$discharge_cfs)
is_peak(dec20$dq[1893:1894])

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
df <- data.frame(Date = dec20$datetime, Q = dec20$discharge_cfs)
df <- ts.format(df)
head(df)

df <- data.frame(Date = flow$dt %>% as.POSIXct, Q = flow$discharge_cfs)
head(df)
ann.cv(df)
baseflows(df)
baseflows(df, ts = 'annual')
Colwells(df)
CTF(df)
daily.cv(df)
flood.length.max(df, 40000)
high.spell.lengths(df)
high.spells(df, threshold = 26000)
partial.series(df, series = T)
seasonality(df, monthly = T)

ps <- partial.series(df[flow$yr > 2000, ], series = T)
psr <- ps$p.series %>% as.data.table
psr <- psr[order(Date), ]
psr
psd <- psr$Date %>% as_datetime
psx <- psd[-1] - psd[-length(psd)]
plot(psr$Q[-1], sort(psx))

length(psx) / sum(psx) %>% as.numeric * 365

dif <- sum(psx)
