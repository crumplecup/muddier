library(data.table)
library(magrittr)


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

rownames(ar) <- ids
colnames(ar) <- index

# subset by facies
df <- ar[ct$facies == 'DF', ]
ff <- ar[ct$facies == 'FF', ]
fg <- ar[ct$facies == 'FG', ]

nrow(df)
nrow(ff)
nrow(fg)


# estimate deposit age by subtracting inherited age fit via convolution

index <- char_pmfs %>% rownames %>% as.numeric

begin <- Sys.time()
da_df <- apply(df, 1, function (x) convo(x, rev(df_ph_pmf), index))
rownames(da_df) <- index
colnames(da_df) <- rownames(df)
end <- Sys.time()
end - begin
# nrow(df) / as.numeric(end - begin)

begin <- Sys.time()
da_ff <- apply(ff, 1, function (x) convo(x, rev(ff_ph_pmf), index))
rownames(da_ff) <- index
colnames(da_ff) <- rownames(ff)
end <- Sys.time()
end - begin
# nrow(ff) / as.numeric(end - begin)

begin <- Sys.time()
da_fg <- apply(fg, 1, function (x) convo(x, rev(fg_ph_pmf), index))
rownames(da_fg) <- index
colnames(da_fg) <- rownames(fg)
end <- Sys.time()
end - begin
# nrow(fg) / as.numeric(end - begin)


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

# remove GRC and CCF samples

rems <- c(grep('GRC', colnames(dadf)), grep('CC', colnames(dadf)))
dadfs <- dadf[,-rems]
rems <- c(grep('GRC', colnames(daff)), grep('CC', colnames(daff)), grep('T8', colnames(daff)))
daffs <- daff[,-rems]
rems <- c(grep('GRC', colnames(dafg)), grep('CC', colnames(dafg)), grep('T8', colnames(dafg)))
dafgs <- dafg[,-rems]

# debris-flow deposit delivery probability at nodes with charcoal samples

creeks_radio$sid[21] <- 'BC-18'
dfids <- colnames(dadfs)
dfdp <- 0
for (i in seq_along(dfids)) {
  dfdp[i] <- creeks$DebrisFlow[creeks$NODE_ID == creeks_radio$node_ids[
    creeks_radio$sid == dfids[i]
  ]]
}

creeks_radio$sid[22] <- 'BC-20'
ffids <- colnames(daffs)
ffdp <- 0
for (i in seq_along(ffids)) {
  ffdp[i] <- creeks$DebrisFlow[creeks$NODE_ID == creeks_radio$node_ids[
    creeks_radio$sid == ffids[i]
    ]]
}

creeks_radio$sid[259] <- 'DFK_184b'
creeks_radio$sid[239] <- 'DFK_140n'
creeks_radio$sid[75] <- 'LK140'
creeks_radio$sid[221] <- 'DFk_42'

fgids <- colnames(dafgs)
fgdp <- 0
for (i in seq_along(fgids)) {
  fgdp[i] <- creeks$DebrisFlow[creeks$NODE_ID == creeks_radio$node_ids[
    creeks_radio$sid == fgids[i]
    ]]
}

# weight delivery probabilities by optimization function df_wt

wdp <- 0
for (i in seq_along(dfdp)) {
  k <- 1
  if (dfdp[i] <= df_wt$dp[k]) {
    wdp[i] <- dfdp[i] * df_wt$wt[k]
  }
  if (dfdp[i] > df_wt$dp[k]) {
    while(dfdp[i] > df_wt$dp[k+1] & k < length(dfdp)) {
      k <- k + 1
    }
    wdp[i] <- dfdp[i] * df_wt$wt[k]
  }
}

wdp <- wdp / sum(wdp)

weighter <- function(pmfs, wts) {
  wpmfs <- array(0, dim(pmfs))
  for (i in seq_along(wts)) {
    vec <- 0
    for (j in seq_along(wts)) {
      vec[j] <- abs(wts[j] - wts[i])
    }
    vec <- vec / max(vec)
    ar <- mapply(function(x, y, z) z[,x] * (1 - y), 1:ncol(pmfs), vec, MoreArgs = list(z = pmfs))
    sumar <- apply(ar, 1, sum)
    wpmfs[ , i] <- sumar / sum(sumar)
  }
  excs <- apply(wpmfs, 2, to_exceed)
  hzd <- array(0, dim(pmfs))
  for (i in 1:nrow(wpmfs)) {
    if (i == 1) {
      hzd[i,] <- wpmfs[i,]
    }
    if (i > 1) {
      hzd[i,] <- wpmfs[i,] / excs[i-1,]
    }
  }
  # hzd[hzd > 1] <- 1
  # hzd[hzd < 0] <- 0
  return(hzd)
}
# for each df pmf
# get the absolute wdp distance from all other pmfs
# weight df pmfs by wdp distance
# sum pmfs and renormalize

whzd <- weighter(dadfs, wdp)

dfmn <- apply(dadfs, 2, function(x) weighted.mean(as.numeric(rownames(dadfs)), x))
ffmn <- apply(dadfs, 2, function(x) weighted.mean(as.numeric(rownames(daffs)), x))
fgmn <- apply(dadfs, 2, function(x) weighted.mean(as.numeric(rownames(dafgs)), x))

plot(dfmn, xlab = 'observation', ylab = 'transit time (yrs)')
lines(cumsum(rexp(ncol(dadfs), rt)), col = 'slateblue')
lines(cumsum(rexp(ncol(dadfs), rt)), col = 'goldenrod')
lines(cumsum(rexp(ncol(dadfs), rt)), col = 'forestgreen')
lines(cumsum(sort(rexp(ncol(dadfs), rt))), col = 'slateblue')
lines(cumsum(sort(rexp(ncol(dadfs), rt))), col = 'goldenrod')
lines(cumsum(sort(rexp(ncol(dadfs), rt))), col = 'forestgreen')


# synthesize inventory of arrivals and removals

rt <- ncol(dadfs) / max(c(dfmn, ffmn, fgmn))
site_n <- ncol(dadfs) + ncol(daffs) + ncol(dafgs)

rtin <- rt
ins <- rexp(ncol(dadfs), rt)
outs <- rexp(ncol(dadfs), rt/10)
plot(sort(cumsum(outs)))
points(sort(cumsum(ins)), col = 'slateblue')

uniq_sam <- function(pool, ls)  {
  s <- sample(pool, 1)
  if (s %in% ls) {
    uniq_sam(pool, ls)
  }
  return(s)
}

recorder <- function(n, ri, ro)  {
  im <- 0
  k <- 1
  rec <- 0
  om <- 0

  while (length(rec) < n)  {
    om <- om + rexp(1, ro)
    while (im < om & length(rec) <= n + 1) {
      im <- im + rexp(1, ri)
      rec[k] <- im
      k <- k + 1
    }
    m <- length(rec[rec < om])
    if (m >= 1 ) {
      s <- round(runif(1, 1, m))
      rec <- rec[-s]
      k <- k - 1
    }
  }
  rec <- rec[1:n]
  return(sort(max(rec) - rec))
}

plot(recorder(ncol(dadfs), rt*10, rt*9.58))
recorder(ncol(dadfs), rt*10, rt*9.58) %>% length



rec <- recorder(1000, rt*10, rt*9.58)
plot(dfmns, pch = 20, col = 'forestgreen', xlab = 'samples', ylab = 'transit time')
points(sort(max(rec[[3]]) - rec[[3]]))
legend('topleft', legend = c('observed', 'synthetic'), fill = c('forestgreen', 'black'))

# compare df cdf to synthetic cdfs
dadfcdf <- apply(dadfs, 1, function(x) sum(x) / length(x)) %>% to_cdf
ogpmfs <- char_pmfs %>% as.matrix
ogpmfs <- ogpmfs[ , colnames(ogpmfs) %in% colnames(dadfs)]
dfcdf <- apply(ogpmfs, 1, function(x) sum(x) / length(x)) %>% rev %>% to_cdf
plot(dfcdf)
discdf <- apply(ogpmfs, 2, function(x) weighted.mean(as.numeric(rownames(char_pmfs)), x)) %>% emp_cdf


boot_rec <- function(n, ri, ro, boot)  {
  ar <- array(0, c(n, boot))
  for (i in 1:boot) {
    ar[, i] <- recorder(n, ri, ro)
  }
  vec <- apply(ar, 1, median)
  return(vec)
}

library(gtools)
vec <- seq(.001, .1, .001)
com <- combinations(length(vec), 2, vec)



cdf_gof <- function(ob1, ob2)  {
  vals <- sort(unique(c(ob1, ob2)))
  c1 <- 0
  c2 <- 0
  for (i in seq_along(vals)) {
    c1[i] <- length(ob1[ob1 <= vals[i]]) / length(ob1)
    c2[i] <- length(ob2[ob2 <= vals[i]]) / length(ob2)
  }
  ks <- max(abs(c1 - c2))  # kolmogorov-smirnov test
  kp1 <- max(c1 - c2)
  kp2 <- max(c2 - c1)
  kp <- kp1 + kp2  # kuiper statistic
  return(data.frame(ks = ks, kp = kp))
}

bit <- cdf_gof(brec[,1], emp_cdf(dfmns)[,1])

vec <- seq(.001, .2, .001)
com <- combinations(length(vec), 2, vec)

# for each combo of rates, bootstrap and get gof stats
gofs <- array(0, c(nrow(com), 2))
sz <- ncol(ogpmfs)
bts <- 100
bcdf <- emp_cdf(dfmns)[,1] %>% as.numeric
dfmns <- apply(ogpmfs, 2, function(x) weighted.mean(rev(index), x))
begin <- Sys.time()
for (i in 1:nrow(com)) {
  b <- boot_rec(sz, com[i,2], com[i,1], bts)
  gofs[i,] <- cdf_gof(b, dfmns) %>% as.numeric
  print(i)
}
end <- Sys.time()
end - begin

scatter3D(com[,1], com[,2], gofs[,1])
scatter3D(com[,1], com[,2], gofs[,2])
ksmin <- com[gofs[,1] == min(gofs[,1]),]
kpmin <- com[gofs[,2] == min(gofs[,2]),]
mrg <- apply(gofs, 1, sum)
mrgmin <- com[mrg == min(mrg),]
scatter3D(com[,1], com[,2], mrg)

png('fit_io.png', height = 14, width = 17, units = 'cm', res = 300)
plot(emp_cdf(dfmns),
     xlab = 'transit time', ylab = 'cdf')
for (i in 1:nrow(ksmin)) {
  lines(emp_cdf(boot_rec(sz, ksmin[i,2], ksmin[i,1], bts)),
        lwd = 2, col = get_palette('crimson', .5))
}
for (i in 1:nrow(kpmin)) {
  lines(emp_cdf(boot_rec(sz, kpmin[i,2], kpmin[i,1], bts)),
        lwd = 2, col = get_palette('ocean', .5))
}
legend('bottomright', legend = c('observed', 'K-S fit', 'Kuiper fit'),
       fill = get_palette(c('charcoal', 'crimson', 'ocean'), .7))
dev.off()


bfs_io <- data.frame(ri = com[,2],
                     ro = com[,1],
                     ks = gofs[,1],
                     kp = gofs[,2])
save(bfs_io, file = 'bfs_io.rds')

boot_exp <- function(sz, rt, bt) {
  ar <- array(0, c(sz, bt))
  for (i in 1:bt) {
    ar[, i] <- cumsum(rexp(sz, rt))
  }
  vec <- apply(ar, 1, median)
  return(vec)
}

png('eg_cdf1.png', height = 14, width = 17, units = 'cm', res = 300)
plot(index, dfcdf, xlim = c(0,10000), xlab = 'transit time', ylab = 'cdf',
     type = 'l', lwd = 2, col = get_palette('charcoal', .5))
lines(index, dadfcdf, lwd = 2, col = get_palette('ocean', .7))
points(emp_cdf(dfmns), pch = 20, col = get_palette('charcoal', .5))
lines(emp_cdf(cumsum(rexp(sz, ncol(ogpmfs) / max(dfmns)))), lwd = 2, col = get_palette('crimson', .7))
legend('bottomright', legend = c('charcoal', 'corrected', 'synthetic'),
       fill = get_palette(c('charcoal', 'ocean', 'crimson'), .7))
dev.off()

rt <- ncol(ogpmfs) / max(dfmns)

png('eg_cdf2.png', height = 14, width = 17, units = 'cm', res = 300)
plot(index, dfcdf, xlim = c(0,10000), xlab = 'transit time', ylab = 'cdf',
     type = 'l', lwd = 2, col = get_palette('charcoal', .5))
lines(index, dadfcdf, lwd = 2, col = get_palette('ocean', .7))
points(emp_cdf(dfmns), pch = 20, col = get_palette('charcoal', .5))
lines(emp_cdf(boot_exp(sz, rt, bts)), lwd = 2, col = get_palette('crimson', .7))
for (i in 1:5) {
  lines(emp_cdf(boot_rec(sz, rt*(i+1), rt*i, bts)), lwd = 2, col = get_palette('crimson', .7))
}
legend('bottomright', legend = c('charcoal', 'corrected', 'synthetic'),
       fill = get_palette(c('charcoal', 'ocean', 'crimson'), .7))
text(8300, .87, 'rate out = 0')
text(6200, .87, '1')
text(4400, .87, '2')
text(3300, .87, '3')
text(2500, .87, '4')
text(2000, .87, '5')
dev.off()


plot(emp_cdf(rexp(sz, rt)),
     xlab = 'interarrival time', ylab = 'cdf')

png('bfs_io.png', height = 20, width = 25, units = 'cm', res = 300)
scatter3D(bfs_io[,2], bfs_io[,1], bfs_io[,3], ticktype = 'detailed',
          xlab = 'output rate', ylab = 'input rate', zlab = 'K-S stat')
dev.off()


brec <- function(n, ri, ro, ia, it = 100) {
  ar <- array(0, c(n, it))
  for (i in 1:it) {
    ar[ , i] <- reco(n, ri, ro, ia) %>% sort
  }
  vec <- apply(ar, 1, median)
  return(vec)
}

plot(emp_cdf(brec(sz, ksmin[1,2], ksmin[1,1], ia_ar)))

reco <- function(n, ri, ro, ia)  {
  im <- 0
  k <- 1
  rec <- 0
  om <- 0

  while (length(rec) < n)  {
    om <- om + rexp(1, ro)
    while (im < om & length(rec) <= n + 1) {
      im <- im + rexp(1, ri)
      rec[k] <- im
      k <- k + 1
    }
    m <- length(rec[rec < om])
    if (m >= 1 ) {
      s <- round(runif(1, 1, m))
      rec <- rec[-s]
      k <- k - 1
    }
  }
  rec <- rec[1:n]
  rec <- sort(max(rec) - rec)
  rc <- 0
  for (i in seq_along(rec)) {
    rc[i] <- rec[i] + sample(ia, 1)
  }
  return(rc - min(rc))
}

plot(emp_cdf(reco(sz, .55, .33, ia_ar)))
reco(sz, .55, .33, ia_ar) %>% sort
# for debris flow deposits with multiple samples
# convolve the youngest age from all ages
# take the estimated mean arrival time of inherited charcoal ages
# estimate lambda based on arrivals in time frame
index <- char_pmfs %>% rownames %>% as.numeric
df_ranks <- rank_list('DF')
df_cl <- convo_list(df_ranks)
df_trc <- lapply(df_cl, function(a) trunc_list(a, sort(index)))
df_trc <- lapply(df_trc, rack)
df_pmf <- df_trc[[1]]

for (i in 2:7) {
  df_pmf <- cbind(df_pmf, df_trc[[i]])
}

ia_ar <- apply(df_pmf, 2, function(x) weighted.mean(rev(index), x))
ia_rt <- length(ia_ar) / max(ia_ar)
ia_cdf <- emp_cdf(ia_ar)[,2]

ia_gofs <- array(0, c(nrow(com), 2))

begin <- Sys.time()
for (i in 1:nrow(com)) {
  b <- brec(sz, com[i,2], com[i,1], ia_ar, bts)
  ia_gofs[i,] <- cdf_gof(b, dfmns) %>% as.numeric
  print(i)
}
end <- Sys.time()
end - begin

scatter3D(com[,1], com[,2], ia_gofs[,1])
scatter3D(com[,1], com[,2], ia_gofs[,2])
ksmin <- com[ia_gofs[,1] == min(ia_gofs[,1]),]
kpmin <- com[ia_gofs[,2] == min(ia_gofs[,2]),]

bfs_ia <- data.frame(ri = com[,2],
                     ro = com[,1],
                     ks = ia_gofs[,1],
                     kp = ia_gofs[,2])
save(bfs_ia, file = 'bfs_ia.rds')

png('bfs_ia.png', height = 20, width = 25, units = 'cm', res = 300)
scatter3D(bfs_ia[,2], bfs_ia[,1], bfs_ia[,3], ticktype = 'detailed', pch = 20,
          xlab = 'output rate', ylab = 'input rate', zlab = 'K-S stat')
dev.off()

bfs_ia[bfs_ia$ks == min(bfs_ia$ks),]
bfs_ia[bfs_ia$kp == min(bfs_ia$kp),]

bts <- 1000
png('fit_ia.png', height = 14, width = 17, units = 'cm', res = 300)
plot(emp_cdf(dfmns),
     xlab = 'transit time', ylab = 'cdf')
for (i in 1:nrow(ksmin)) {
  lines(emp_cdf(brec(sz, ksmin[1,2], ksmin[1,1], ia_ar, bts)),
        lwd = 2, col = get_palette('crimson', .33))
}
for (i in 1:nrow(kpmin)) {
  lines(emp_cdf(brec(sz, kpmin[i,2], kpmin[i,1], ia_ar, bts)),
        lwd = 2, col = get_palette('ocean', .33))
}
legend('bottomright', legend = c('observed', 'K-S fit', 'Kuiper fit'),
       fill = get_palette(c('charcoal', 'crimson', 'ocean'), .7))
dev.off()

meander <- function(n, ri, ro, ia, scl, mr)  {
  rec <- 0
  plc <- 0
  ctr <- rnorm(1, sd = scl)
  im <- 0
  om <- 0
  k <- 1

  while (length(rec) < n)  {
    om <- om + rexp(1, ro)
    # print(paste0('Removal on year ', round(om, 2)))
    while (im < om & length(rec) <= n + 1) {
      im <- im + rexp(1, ri)
      rec <- c(rec, im)
      # print('Arrival year ')
      # print(round(im, 2))
      if (runif(1) < mr) ctr <- rnorm(1, sd = scl)
      plc <- c(plc, ctr)
      # print(paste0('k = ', k))
      # print('Places ')
      # print(plc)
      # k <- k + 1
    }
    m <- length(rec[rec < om])
    if (length(m) > 1 ) {
      # print('Valid Places ')
      # print(plc[rec < om])
      d <- abs(ctr - plc[rec < om])  # absolute distance from center
      # print('Absolute distances ')
      # print(d)
      d <- 1 - d / max(d)  # inverse distance proportion
      # print('Inverse distance proportion')
      # print(d)
      d <- (1 / length(d)) * d  # multiply by srs probability
      # print('weighted probability')
      # print(d)
      d <- d / sum(d) # renormalize
      # print('normalized pdf')
      # print(d)
      d <- cumsum(d)  # convert to cdf
      # print('cdf')
      # print(d)
      r <- runif(1)
      rmv <- abs(r - plc[rec < om])
      rmv <- sum(as.numeric(rmv == min(rmv) * 1:length(rmv)))
      # print(paste0("Removing ", rmv))
      rec <- rec[-rmv]
      plc <- plc[-rmv]
      k <- k - 1
    }
  }
  rec <- rec[1:n]
  rec <- sort(max(rec) - rec)
  rc <- 0
  for (i in seq_along(rec)) {
    rc[i] <- rec[i] + sample(ia, 1)
  }
  return(rc - min(rc))
}

plot(emp_cdf(meander(sz, ksmin[1,2], ksmin[1,1], ia_ar, 1, .3)))

accumulate <- function(n, ri, ro, ia, og, it = 1000, scl = 1, mr = 1) {
  ar <- array(0, c(n, it))
  for (i in 1:it) {
    ar[ , i] <- meander(n, ri, ro, ia, scl, mr) %>% sort
  }
  med <- apply(ar, 1, median)
  cdf <- emp_cdf(med)
  gofs <- apply(ar, 2, function(x) cdf_gof(x, cdf)) %>% rackl %>% matrix(ncol = 2)
  gof <- cdf_gof(og, cdf)
  ks <- length(gofs[gofs[,1] <= gof[1], 1]) / nrow(gofs)
  kp <- length(gofs[gofs[,2] <= gof[2], 2]) / nrow(gofs)

  return(list(med, data.frame(ks = ks, kp = kp)))
}

accumulate(sz, kpmin[1,2], kpmin[1,1], ia_ar, dfmns, 100)
plot(emp_cdf(accumulate(sz, .55, .33, ia_ar, dfmns, 100)[[1]]))


# redoing the brute force search after fixing some bugs
bts <- 200
ia_gofs <- array(0, c(nrow(com), 2))

begin <- Sys.time()
for (i in 1:nrow(com)) {
  b <- brec(sz, com[i,2], com[i,1], ia_ar, bts)
  ia_gofs[i,] <- cdf_gof(b, dfmns) %>% as.numeric
  print(i)
}
end <- Sys.time()
end - begin

me_gofs <- array(0, c(nrow(com), 4))

begin <- Sys.time()
ln <- nrow(com)
for (i in 1:ln) {
  ac <- accumulate(sz, com[i,2], com[i,1], ia_ar, dfmns, 100)
  acob <- ac[[1]]
  acgof <- ac[[2]] %>% unlist %>% as.numeric
  me_gofs[i, 1:2] <- cdf_gof(acob, dfmns) %>% as.numeric
  me_gofs[i, 3] <- acgof[1]
  me_gofs[i, 4] <- acgof[2]
  tracker(i, ln, begin)
}
end <- Sys.time()
end - begin

scatter3D(com[,1], com[,2], me_gofs[,1])
scatter3D(com[,1], com[,2], me_gofs[,2])
scatter3D(com[,1], com[,2], me_gofs[,3])
scatter3D(com[,1], com[,2], me_gofs[,4])
ksmin <- com[me_gofs[,1] == min(me_gofs[,1]),]
kpmin <- com[me_gofs[,2] == min(me_gofs[,2]),]

bfs_me <- data.frame(ri = com[,2],
                     ro = com[,1],
                     ks = me_gofs[,1],
                     kp = me_gofs[,2])
save(bfs_me, file = 'bfs_me.rds')

ksmin <- bfs_ia[bfs_ia[,3] == min(bfs_ia[,3]),]
kpmin <- bfs_ia[bfs_ia[,4] == min(bfs_ia[,4]),]
bfs_me[bfs_me[,3] == min(bfs_me[,3]),]
bfs_me[bfs_me[,4] == min(bfs_me[,4]),]

png('ks_kp_mins_ia.png', height = 15, width = 17, units = 'cm', res = 300)
plot(ksmin$ro, ksmin$ri, pch = 20, col = get_palette('crimson', .7),
     xlim = c(min(ksmin$ro), max(kpmin$ro)),
     ylim = c(min(ksmin$ri), max(kpmin$ri)),
     xlab = 'output rate', ylab = 'input rate')
points(kpmin$ro, kpmin$ri, pch = 20, col = get_palette('ocean', .7))
legend('bottomright', legend = c('K-S', 'Kuiper'),
       fill = get_palette(c('crimson', 'ocean'), .7))
dev.off()


png('fit_me.png', height = 14, width = 17, units = 'cm', res = 300)
plot(emp_cdf(dfmns),
     xlab = 'transit time', ylab = 'cdf')
for (i in 1:nrow(ksmin)) {
  lines(emp_cdf(accumulate(sz, ksmin[i,2], ksmin[i,1], ia_ar, dfmns, bts)[[1]]),
      lwd = 2, col = get_palette('crimson', .5))
}
lines(emp_cdf(accumulate(sz, kpmin[2], kpmin[1], ia_ar, dfmns, bts)[[1]]),
        lwd = 3, col = get_palette('ocean', .5))

legend('bottomright', legend = c('observed', 'K-S fit', 'Kuiper fit'),
       fill = get_palette(c('charcoal', 'crimson', 'ocean'), .7))
dev.off()

png('bfs_me.png', height = 20, width = 25, units = 'cm', res = 300)
scatter3D(bfs_me[,2], bfs_me[,1], bfs_me[,3], ticktype = 'detailed', pch = 20,
          xlab = 'output rate', ylab = 'input rate', zlab = 'K-S stat')
dev.off()




# older code

scl <- 1
brec <- boot_rec(ncol(dadfs), rt*(scl+1), rt*scl, 20) %>% emp_cdf
brec[,2] %>% ks_test(discdf[,2])
plot(index, dfcdf, xlim = c(0, 10000), type = 'l', lwd = 3,
     col = get_palette('ocean', .7), ylim = c(0,1))
for (i in 1:10) {
  lines(emp_cdf(boot_rec(ncol(dadfs), rt*(i+1), rt*i, 20)), col = get_palette('charcoal', .5))
}
rec <- recorder(1000, rt*10, rt*9.59)
plot(dfmns, pch = 20, col = 'forestgreen', xlab = 'samples', ylab = 'transit time')
points(sort(max(rec[[3]]) - rec[[3]]))
legend('topleft', legend = c('observed', 'synthetic'), fill = c('forestgreen', 'black'))



for (i in 1:10) {
  lines(recorder(ncol(dadfs), rt, rt/i, 10000))
}

ins <- cumsum(rexp(ncol(dadfs), rt))
outs <- cumsum(rexp(ncol(dadfs), rt/2))
k <- 1
rec <- ins
while (outs[k] < t) {
  m <- red[rec <= outs[k]]
  if (length(m) >= 1 ) {
    m <- sample(m, 1)
  }
  rec <- rec[-match(m, rec)]
  k <- k + 1
}





rt <- ncol(dadfs) / max(c(dfmn, ffmn, fgmn))
prds <- rexp(20, .2)
prda <- cumsum(prds)
prda %>% to_hzd
char_hzd <- unlist(charcoal$mn) %>% to_hzd
plot(char_hzd[,4], char_hzd[,1])

dfpmf <- apply(df, 2, function(x) sum(x) / length(x)) %>% rev
dfhzd <- dfpmf / (1 - to_cdf(dfpmf))
dfhzd[dfhzd < 0] <- 0
dfhzd[dfhzd > 1] <- 1

plot(index+1, dfhzd, log = 'xy', type = 'l', lwd = 2)
(dadfs[,40] / (1 - (dadfs[,40] %>% to_cdf))) %>% plot(index, ., log = 'xy')
lmda <- dfhzd[2:100] %>% mean

prds <- rexp(10, lmda)
prda <- cumsum(prds)
prd_cdf <- prda %>% emp_cdf
prd_pmf <- prd_cdf[,2] %>% to_pmf
prd_hzd <- to_hzd(prda)
plot(prd_cdf[,1], prd_hzd)
lines(prd_cdf[,1], prd_cdf[,2], col = 'forestgreen')
abline(h = lmda, lwd = 2, col = 'slateblue', lty = 2)
data.frame(age = prd_cdf[,1], cdf = prd_cdf[,2], hzd = prd_hzd)

ob <- round(runif(50, 1, 6))
plot(emp_cdf(ob)[,1], to_hzd(ob)[,1])
ob %>% to_hzd

to_hzd <- function(obs) {
  cdf <- emp_cdf(obs)
  pmf <- to_pmf(cdf[,2])
  hzd <- 0
  for (i in seq_along(pmf)) {
    if (i == 1) {
      hzd[i] <- pmf[i]
    }
    if (i > 1) {
      hzd[i] <- pmf[i] / (1 - cdf[i-1, 2])
    }
  }
  df <- data.frame(hzd = hzd, pmf = pmf, cdf = cdf[,2], val = cdf[,1])
  return(df)
}


to_pmf5 <- function(obs)  {
  hi <- ceiling(max(obs))
  top <- 5
  while (top < hi) {
    top <- top + 5
  }
  cdf <- emp_cdf(obs)
  cdfs <- vector((top / 5) + 1, mode = 'numeric')
  k <- 0
  for (i in 1:nrow(cdf)) {
    while (cdf[i,1] > k*5) {
      k <- k + 1
      if (k > 1) cdfs[k] <- cdfs[k-1]
    }
    cdfs[k] <- cdf[i,2]
  }
  cdfs[k+1] <- cdfs[k]
  pmf <- to_pmf(cdfs)
  hzd <- vector(length(pmf), mode = 'numeric')
  for (i in seq_along(pmf)) {
    if (i == 1) {
      hzd[i] <- pmf[i]
    }
    if (i > 1) {
      hzd[i] <- pmf[i] / (1 - cdfs[i-1])
    }
  }
  return(hzd)
}
to_pmf5(prda) %>% plot(log = 'y')

studs <- 500
yr1 <- 15
yr2 <- 23
yr2 / (studs - yr1) == yr_hzr[2]
yr_pmf <- c(yr1/studs, yr2/studs, 0)
yr_cdf <- c(0, cumsum(yr_pmf))[-4]
yr_hzr <- yr_pmf / (1 - yr_cdf)

lmda <- whzd[2:100,1] %>% mean
dn <- ncol(dadfs)
prds <- rexp(dn, lmda) %>% emp_cdf
prds_pmf <- prds[,2] %>% to_pmf
prds_hzd <- prds_pmf / (1 - prds[,2])
plot(prds[,1], prds_hzd, log = 'y')

plot(as.numeric(rownames(dadfs)), whzd[,1], log = 'y', type = 'l', lwd = 2, xlim = c(0,10000))
lines(as.numeric(rownames(dadfs)), whzd[,10],
      lwd = 2, col = 'slateblue')
lines(as.numeric(rownames(dadfs)), whzd[,60],
      lwd = 2, col = 'forestgreen')


dfmns <- apply(dadfs, 2, function(x) weighted.mean(as.numeric(rownames(dadfs)), x))
ffmns <- apply(daffs, 2, function(x) weighted.mean(as.numeric(rownames(daffs)), x))
fgmns <- apply(dafgs, 2, function(x) weighted.mean(as.numeric(rownames(dafgs)), x))

dfnids <- lapply(dfids, function(x) creeks_radio$node_ids[creeks_radio$sid %in% x]) %>% unlist

ffids <- colnames(daffs)
ffnids <- lapply(ffids, function(x) creeks_radio$node_ids[creeks_radio$sample_id %in% x]) %>% unlist
ffp <- lapply(ffnids, function(x) creeks$DebrisFlow[creeks$NODE_ID %in% x]) %>% unlist

fgids <- colnames(dafgs)
fgnids <- lapply(fgids, function(x) creeks_radio$node_ids[creeks_radio$sample_id %in% x]) %>% unlist
fgp <- lapply(fgnids, function(x) creeks$DebrisFlow[creeks$NODE_ID %in% x]) %>% unlist

charcrks <- creeks_radio
charcrks$dp <- -1
for (i in 1:nrow(charcrks)) {
  charcrks$dp[i] <- creeks$DebrisFlow[creeks$NODE_ID == charcrks$node_ids[i]]
}

crksdf <- charcrks[charcrks$sample_id == dfids[2],]

dfids %>% unique %>% length

for (i in seq_along(dfids)) {
  crksdf <- rbind(crksdf, charcrks[charcrks$node_ids %in% dfnids[i],])
}
crksdf <- crksdf[-1, ]

mnda <- 0
k <- 1
for (i in seq_along(dfids)) {
  sid <- creeks_radio$sample_id[creeks_radio$node_ids %in% dfnids[i]]
  for (j in 1:length(sid)) {
    pmf <- dadf[ , colnames(dadf) %in% sid[j]]
    mnda[k] <- weighted.mean(pmf, index)
    k <- k + 1
  }
}

crksdf$mnda <- apply(dadfs, 2, function(x) weighted.mean(x, index))
plot(crks$dfp, crks$DebrisFlow)

# the rate at which the interarrival period increases represents
# the removal rate of the stream
# we can represent the removal rate using the hazard functin of deposit ages
# the hazard function is the pmf divided by the exceedance function

dadf_pmf <- apply(dadf, 1, function(x) sum(x) / length(x))
dadf_exc <- to_exceed(dadf_pmf)
dadf_hzd <- dadf_pmf / dadf_exc

daff_pmf <- apply(daff, 1, function(x) sum(x) / length(x))
daff_exc <- to_exceed(daff_pmf)
daff_hzd <- daff_pmf / daff_exc

dafg_pmf <- apply(dafg, 1, function(x) sum(x) / length(x))
dafg_exc <- to_exceed(dafg_pmf)
dafg_hzd <- dafg_pmf / dafg_exc


index <- char_pmfs %>% rownames %>% as.numeric %>% sort
options(scipen = 5)
setwd('/home/crumplecup/work/')

pal <- get_palette(c('crimson', 'charcoal', 'ocean'), .8)
png('da_hzd.png', height = 17, width = 29, units = 'cm', res = 300)
plot(index, dafg_hzd, log = 'xy', type = 'l', lwd = 2, col = pal[3],
     xlab = 'Deposit Age', ylab = 'Hazard Function')
lines(index, daff_hzd, lty = 2, lwd = 2, col = pal[2])
lines(index, dadf_hzd, pch = 20, lwd = 2, col = pal[1])
legend('topleft', legend = c('Debris Flows', 'Fluvial Fines', 'Fluvial Gravels'),
       fill = get_palette(c('crimson', 'charcoal', 'ocean'), .8))
dev.off()
