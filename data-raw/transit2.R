library(data.table)
library(magrittr)
library(parallel)
library(plot3D)
options(mc.cores=detectCores())
set.seed(10101)


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

write.csv(dfmn, 'dfmn.csv')
write.csv(ffmn, 'ffmn.csv')
write.csv(fgmn, 'fgmn.csv')


options(scipen = 5)
pal = get_palette(c('forest', 'coral', 'ocean'))
# png('ia_da.png', height = 12, width = 17, units = 'cm', res = 300)
plot(dfmn[-1], itmn, log = 'x', pch = 20, col = get_palette('slate'),
     xlab = 'Deposit Age', ylab = 'Interarrival Time')
points(fgmn[-1], fgitmn, col = pal[3], pch = 20)
points(ffmn[-1], ffitmn, col = pal[2], pch = 20)
points(dfmn[-1], itmn, col = pal[1], pch = 20)
legend('topleft', legend = c('debris flows', 'fluvial fines', 'fluvial gravels'),
       fill = get_palette(c('forest', 'coral', 'ocean'), .7))
# dev.off()

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


# debris-flow deposit delivery probability at nodes with charcoal samples

creeks_radio$sid[21] <- 'BC-18'
dfdp <- 0
for (i in seq_along(dfids)) {
  dfdp[i] <- creeks$DebrisFlow[creeks$NODE_ID == creeks_radio$node_ids[
    creeks_radio$sid == dfids[i]
  ]]
}

creeks_radio$sid[22] <- 'BC-20'
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


# for each df pmf
# get the absolute wdp distance from all other pmfs
# weight df pmfs by wdp distance
# sum pmfs and renormalize

dfmns <- apply(dfs, 1, function(x) weighted.mean(as.numeric(colnames(dfs)), x))
ffmns <- apply(ffs, 1, function(x) weighted.mean(as.numeric(colnames(ffs)), x))
fgmns <- apply(fgs, 1, function(x) weighted.mean(as.numeric(colnames(fgs)), x))

setwd('/home/crumplecup/work/')
write.csv(data.frame(age = dfmns), 'dfmns.csv')
write.csv(data.frame(age = ffmns), 'ffmns.csv')
write.csv(data.frame(age = fgmns), 'fgmns.csv')


# synthesize inventory of arrivals and removals
library(gtools)
vec <- seq(.001, .25, .001)
com <- combinations(length(vec), 2, vec)
idm <- matrix(c(vec,vec), ncol = 2)
m <- rbind(com, rbind(idm, com[ , c(2,1)]))

uniq_sam <- function(pool, ls)  {
  s <- sample(pool, 1)
  if (s %in% ls) {
    uniq_sam(pool, ls)
  }
  return(s)
}


boot_exp <- function(sz, rt, bt) {
  ar <- array(0, c(sz, bt))
  for (i in 1:bt) {
    ar[, i] <- cumsum(rexp(sz, rt))
  }
  vec <- apply(ar, 1, median)
  return(vec)
}



# for each combo of rates, bootstrap and get gof stats

# for debris flow deposits with multiple samples
# convolve the youngest age from all ages
# take the estimated mean arrival time of inherited charcoal ages
# use the inherited ages as part of the bootstrap
index <- char_pmfs %>% rownames %>% as.numeric
df_ranks <- rank_list('DF')
df_cl <- convo_list(df_ranks)
df_trc <- lapply(df_cl, function(a) trunc_list(a, sort(index)))
df_trc <- lapply(df_trc, rack)
df_pmf <- df_trc[[1]]

for (i in 2:7) {
  df_pmf <- cbind(df_pmf, df_trc[[i]])
}

ff_ranks <- rank_list('FF')
ff_cl <- convo_list(ff_ranks)
ff_trc <- lapply(ff_cl, function(a) trunc_list(a, sort(index)))
ff_trc <- lapply(ff_trc, rack)
ff_pmf <- ff_trc[[1]]

for (i in 2:3) {
  ff_pmf <- cbind(ff_pmf, ff_trc[[i]])
}

fg_ranks <- rank_list('FG')
fg_cl <- convo_list(fg_ranks)
fg_trc <- lapply(fg_cl, function(a) trunc_list(a, sort(index)))
fg_trc <- lapply(fg_trc, rack)
fg_pmf <- fg_trc[[1]]

for (i in 2:4) {
  fg_pmf <- cbind(fg_pmf, fg_trc[[i]])
}

ncol(df_pmf) + ncol(ff_pmf) + ncol(fg_pmf) +
nrow(df) + nrow(ff) + nrow(fg)
ia_pmf <- cbind(df_pmf, ff_pmf, fg_pmf)
dfia <- apply(df_pmf, 2, function(x) weighted.mean(rev(index), x))
ffia <- apply(ff_pmf, 2, function(x) weighted.mean(rev(index), x))
fgia <- apply(fg_pmf, 2, function(x) weighted.mean(rev(index), x))
dfia <- data.frame(ia = dfia)
dfia$type <- 'DF'
ffia <- data.frame(ia = ffia)
ffia$type <- 'FF'
fgia <- data.frame(ia = fgia)
fgia$type <- 'FG'
iat <- rbind(dfia, rbind(ffia, fgia))
write.csv(iat, 'iat.csv')

dfdf <- data.frame(age = dfmn)
dfdf$type <- 'DF'
ffdf <- data.frame(age = ffmn)
ffdf$type <- 'FF'
fgdf <- data.frame(age = fgmn)
fgdf$type <- 'FG'
dep <- rbind(dfdf, rbind(ffdf, fgdf))
write.csv(dep, 'dep.csv')

ia_ar <- apply(ia_pmf, 2, function(x) weighted.mean(rev(index), x))
plot(emp_cdf(ia_ar))

write.csv(ia_ar, file = 'ia_ar.csv')
# resource consuming bootstrap
bts <- 10
sz <- nrow(dfs)
ln <- nrow(m)
ia_gofs <- array(0, c(nrow(com), 2))

boot_wrap <- function(n, ri, ro, ia, it, cdf) {
  mcmapply(function(x,y,a,b,c,d) cdf_gof(accumulate(a, x, y, b, c), as.numeric(d)),
         x = ri, y = ro, MoreArgs = list(a = n, b = ia, c = it, d = cdf))
}

iag <- mcmapply(function(x,y, a,b,c) accumulate(a, x, y, b, c), x = m[,1], y = m[,2],
              MoreArgs = list(a = sz, b = ia_ar, c = bts))

for_wrap <- function(n, ri, ro, ia, it, cdf) {
  ln <- length(ri)
  print(ln)
  gofs <- array(0, dim = c(ln, 2))
  for (i in 1:ln) {
    b <- accumulate(n, ri, ro, ia, it)
    gofs[i,] <- cdf_gof(b, cdf) %>% as.numeric
  }
  gofs
}

library(microbenchmark)

interval <- 19001:19100
microbenchmark(
  wapply = boot_wrap(sz, m[interval,1], m[interval,2], ia_ar, bts, dfmn),
  wfor = for_wrap(sz, m[interval,1], m[interval,2], ia_ar, bts, dfmn)
)

# make sure the inverse cases are doing what you think

bit <- recorder(sz, .01, .01, ia_ar)

bit <- accumulate(sz, m[1000,1], m[1000,2], ia_ar, 100)

bot <- cdf_gof(bit, dfmn)
(bit <- rexp(10, .01))


# too big to bootstrap at once, do in pieces

chnk <- com[com[,1] > .2 & com[,2] > .2, ]
chnk1 <- com[com[,1] <= .2 & com[,2] > .2, ]

bts <- 200 # Time difference of 1.015347 days
gofs <- array(0, c(nrow(chnk), 2))
begin <- Sys.time()
ln <- nrow(chnk)
for (i in 1:ln) {
  b <- accumulate(sz, chnk[i,2], chnk[i,1], ia_ar, bts)
  gofs[i,] <- cdf_gof(b, dfmn) %>% as.numeric
  tracker(i, ln, begin)
}
end <- Sys.time()
end - begin

bfs_ia_chnk <- data.frame(ri = chnk[,2],
                         ro = chnk[,1],
                         ks = gofs[,1],
                         kp = gofs[,2])
save(bfs_ia_chnk, file = 'bfs_ia_chnk.rds')


bts <- 200 #
gofs <- array(0, c(nrow(chnk1), 2))
begin <- Sys.time()
ln <- nrow(chnk1)
for (i in 1:ln) {
  b <- accumulate(sz, chnk1[i,2], chnk1[i,1], ia_ar, bts)
  gofs[i,] <- cdf_gof(b, dfmn) %>% as.numeric
  tracker(i, ln, begin)
}
end <- Sys.time()
end - begin

bfs_ia_chnk1 <- data.frame(ri = chnk1[,2],
                          ro = chnk1[,1],
                          ks = gofs[,1],
                          kp = gofs[,2])
save(bfs_ia_chnk1, file = 'bfs_ia_chnk1.rds')


bts <- 50 # Time difference of 1.015347 days
ia_gofs <- array(0, c(nrow(com), 2))
begin <- Sys.time()
ln <- nrow(com)
for (i in 1:ln) {
  b <- accumulate(sz, com[i,1], com[i,2], ia_ar, bts)
  ia_gofs[i,] <- cdf_gof(b, dfmn) %>% as.numeric
  tracker(i, ln, begin)
}
end <- Sys.time()
end - begin

setwd('/home/crumplecup/work/')
bfs_ia_upr <- data.frame(ri = com[,1],
                     ro = com[,2],
                     ks = ia_gofs[,1],
                     kp = ia_gofs[,2])
save(bfs_ia_upr, file = 'bfs_ia_upr.rds')



bts <- 1000
begin <- Sys.time()
ln <- nrow(idm)
idm_gofs <- array(0, c(nrow(idm), 2))
for (i in 1:ln) {
  b <- accumulate(sz, idm[i,1], idm[i,2], ia_ar, bts)
  idm_gofs[i,] <- cdf_gof(b, dfmn) %>% as.numeric
  tracker(i, ln, begin)
}
end <- Sys.time()
end - begin

scatter3D(idm[,1], idm[,2], idm_gofs[,1])
idm[idm_gofs[,1] == min(idm_gofs[,1]),]

bfs_idm <- data.frame(ri = idm[,1],
                     ro = idm[,2],
                     ks = idm_gofs[,1],
                     kp = idm_gofs[,2])
save(bfs_idm, file = 'bfs_idm.rds')


bfs_io <- rbind(
  bfs_ia, rbind(
    bfs_idm, rbind(
      bfs_ia_upr, rbind(
        bfs_ia_chnk, rbind(
          bfs_ia_chnk1
        )
      )
    )
  )
)


png('bfs_io.png', height = 30, width = 30, units = 'cm', res = 300)
scatter3D(bfs_io[,2], bfs_io[,1], bfs_io[,3], ticktype = 'detailed', pch = 20,
          xlab = 'output rate', ylab = 'input rate', zlab = 'K-S stat',
          phi = 20, theta = 290)
dev.off()

scatter3D(creeks_xyz[,1], creeks_xyz[,2], creeks_xyz[,3],
          theta = 150, pch = 20)

# make CDF of cross-sectional area at nodes with charcoal samples

# subset nodes at bear and knowles with radiocarbon samples

# subset charcoal data for bear and knowles
char <- charcoal$site_id[!charcoal$site_id %in% ct$ids &
                   charcoal$family %in% charcoal$family[ct$N > 1]]

crk_rad <- creeks_radio[creeks_radio$sid %in% char, ]
char[!char %in% crk_rad$sid]
# change site ids to match creeks_radio
# creeks_radio$sid[!creeks_radio$sid %in% crk_rad$sid]
char[char == 'DFK_140m'] <- 'DFK_140M'
char[char == 'DFK_140q'] <- 'DFK_140Q'
char[char == 'DFK_140r'] <- 'DFK_140R'
char[char == 'DFK_140l'] <- 'DFK_140L'
char[char == 'DFK_140p'] <- 'DFK_140P'
char[char == 'DFK_140s'] <- 'DFK_140S'
char[char == 'LK_0.16m'] <- 'LK_0.16M'
char[char == 'LK_0.16l'] <- 'LK_0.16L'
char[char == 'LK_0.42q'] <- 'LK_0.42Q'
char[char == 'LK_118r'] <- 'LK_118R'
char[char == 'LK_151n'] <- 'LK_151N'
char[char == 'LK_145f'] <- 'LK_145F'
char[char == 'LK_145o'] <- 'LK_145O'
char[char == 'LK_145p'] <- 'LK_145P'
char[char == 'LK_145m'] <- 'LK_145M'
char[char == 'LK_187l'] <- 'LK_187L'

# figure out which charcoal sample corresponds to which node in the creeks data
crk <- creeks[creeks$creek_name %in% c('bear', 'knowles'), ]
crk <- crk[crk$NODE_ID %in% creeks_radio$node_ids, ]
sids <- creeks_radio$sid[creeks_radio$creek_name %in% c('bear', 'knowles')]
gids <- 0
bids <- 0
nids <- 0
k <- 1
l <- 1
for (i in 1:nrow(crk)) {
  ids <- creeks_radio$node_ids[creeks_radio$node_ids %in% crk$NODE_ID[i]]
  if (length(ids) == 0) {
    bids[l] <- crk$NODE_ID[i]
    l <- l + 1
  }
  if (length(ids) > 0) {
    fids <- creeks_radio$sid[creeks_radio$node_ids %in% crk$NODE_ID[i]]
    for (j in seq_along(ids)) {
      gids[k] <- fids[j]
      nids[k] <- ids[j]
      k <- k + 1
    }
  }
}

# because some sites sample multiple strata, or single deposits multiple times
# there are 106 unique nodes among 145 samples in bear and knowles

# since we know the input/output rate and optimized del prob for each node
# we can include xsec_area for all nodes in the bear and knowles study area
crks <- creeks[creeks$creek_name %in% c('bear', 'knowles'), ]
setwd('/home/crumplecup/work/')
png('br_kn_xsec_cdf.png', height = 15, width = 17, units = 'cm', res = 300)
plot(crks$xsec_area %>% emp_cdf, pch = 20, col = get_palette('ocean'),
     xlab = 'unit cross-sectional area (m2)', ylab = 'CDF')
lines(crks$xsec_area %>% emp_cdf, lty = 2, col = get_palette('charcoal'))
dev.off()

# make dataframe of study area node coords, slope
crk_coords <- coordinates(crks)
cslp <- crks$slope
bslp <- crks$slope[crks$creek_name == 'bear']
bkm <- crks$ToMouth_km[crks$creek_name == 'bear']
kslp <- crks$slope[crks$creek_name == 'knowles']
kkm <- crks$ToMouth_km[crks$creek_name == 'knowles']

slope_to_elev <- function(slope, dist) {
  elev <- 0
  for (i in seq_along(slope)) {
    if (i == 1) {
      elev[i] <- 0
    }
    if (i > 1) {
    elev[i] <- elev[i-1] + ((dist[i] - dist[i-1]) * slope[i-1])
    }
  }
  elev
}

belev <- slope_to_elev(bslp, bkm)
kelev <- slope_to_elev(kslp, kkm)
elev <- c(belev,kelev)
plot(elev)

creeks_xyz <- data.frame(
  x_data = crk_coords[,1],
  y_data = crk_coords[,2],
  z_data = elev,
  slope = crks$slope
  )

write.csv(creeks_xyz, file = 'creeks_xyz.csv')


# estimate interarrival volume rate from observed xsec areas
# define search space ranging smaller and larger by an order of magnitude
vrt <- nrow(crks) / max(crks$xsec_area)
sq <- seq(.0, 15, .5)
csq <- combinations(length(sq), 2, sq)
idc <- matrix(c(sq,sq), ncol = 2)

crks$lvl <- crks$xsec_area / crks$valley_width
crk$lvl <- crk$xsec_area / crk$valley_width
plot(crks$lvl %>% emp_cdf)



# function to convert MB del probs to optimized del probs

weight_by_delprob <- function(delprob, weight = df_wt) {
  flag <- 0
  k <- 1
  ln <- nrow(weight)
  while (k <= ln & flag == 0) {
    if (weight[k, 2] >= delprob) {
      flag <- k
    }
    if (weight[k, 2] < delprob & k < ln) {
      k <- k + 1
    }
    if (weight[k, 2] < delprob & k == ln) {
      flag <- k
    }
  }

  res <- 0
  if (flag > 0) {
    res <- weight[flag, 1]
  }
  if (flag == 0) {
    res <- NA
  }
  res
}


wt <- lapply(crk$DebrisFlow, function(x) weight_by_delprob(x)) %>% unlist
crk$dp <- crk$DebrisFlow * wt
wt <- lapply(crks$DebrisFlow, function(x) weight_by_delprob(x)) %>% unlist
crks$dp <- crks$DebrisFlow * wt

# global input/output rates
gi <- 0.11
go <- 0.08

# localized weighted rate
# nrow(crk) is the number of sampled nodes in the study areas
# not the total nodes in the study area crks

ri_ave <- gi / nrow(crk)
ri_ttl <- ri_ave * nrow(crks)
ri_dp <- ri_ave * crks$dp
crks$ri <- ri_dp * (ri_ttl / sum(ri_dp))


# output rate is dependent on streampower
# relative streampower based on contr area and slope

crks$pwr <- log(crks$contr_area+1) * log(crks$slope + 1)

pwr_df <- data.frame(CA = log(crks$contr_area),
                     slope = log(crks$slope),
                     creek = crks$creek_name)
pwr_mod <- lm(CA ~ slope + creek, data = pwr_df)
summary(pwr_mod)

pwr_prd <- predict(
  pwr_mod, newdata = data.frame(slope = log(crks$slope), creek = crks$creek_name)
  )

creek_bool <- rep(0, nrow(crks))
creek_bool[crks$creek_name == 'knowles'] <- 1

lca_obs <- log(crks$slope) * pwr_mod$coefficients[2] +
  creek_bool * pwr_mod$coefficients[3]

lca_dif <- log(crks$contr_area) - pwr_prd
difs <- (exp(lca_dif))
# values of difs from zero to one had negative logged difs
# no need to normalize to use as a weight, as negatives are deflated, positive inflated
# streampower coefficient (difs) is normalized to use as a weight

# output weighted by inverse dp and streampower coefficient
ro_ave <- go / nrow(crk)
ro_ttl <- ro_ave * nrow(crks)
ro_wt <- (1 - crks$dp / max(crks$dp)) * difs
crks$ro <- ro_wt * (ro_ttl / sum(ro_wt))
crks$ro[crks$ro <= 0] <- 0.00000001

plot(crks$ri, crks$ro)
points(crks$ri[crks$ri > crks$ro], crks$ro[crks$ri > crks$ro], col = 'red')
plot(crks$ri - crks$ro, pch = 20, col = get_palette('charcoal'))
abline(h = 0, lty = 2)

# function to record volumes
rec_vol <- function(node, n, vi, vo) {
  t <- 0
  age <- 0
  vol <- 0
  k <- 0
  while (t < n) {
    k <- k + 1
    t <- t + rexp(1, node$ri)
    age[k] <- t
    vol[k] <- rexp(1, vi)
  }
  df <- data.frame(
    t = age,
    r = 'input',
    vol = vol,
    lvl = vol / as.numeric(node$valley_width)
    )
  t <- 0
  age <- 0
  vol <- 0
  k <- 0
  while (t < n) {
    k <- k + 1
    t <- t + rexp(1, node$ro)
    age[k] <- t
    vol[k] <- rexp(1, vo)
  }
  df <- rbind(df, data.frame(
    t = age,
    r = 'output',
    vol = vol,
    lvl = vol / as.numeric(node$valley_width)))
  df <- df[order(df$t),]
  vol <- 0
  lvl <- 0
  for (i in 1:nrow(df)) {
    if (df$r[i] == 'input') {
      vol <- vol + df$vol[i]
      lvl <- lvl + df$lvl[i]
    }
    if (df$r[i] == 'output') {
      if (df$vol[i] >= vol) {
        vol <- 0
      }
      if (df$vol[i] < vol) {
        vol <- vol - df$vol[i]
      }
      if (df$lvl[i] >= lvl) {
        lvl <- 0
      }
      if (df$lvl[i] < lvl) {
        lvl <- lvl - df$lvl[i]
      }
    }
  }
  data.frame(vol = vol, lvl = lvl)
}




# function to fit volume rate to observed volumes

fit_volumes <- function(nodes, n, vi, vo, it = 10) {
  vol <- array(0, c(nrow(nodes), it))
  lvl <- array(0, c(nrow(nodes), it))
  for (i in 1:it) {
    for (j in 1:nrow(nodes)) {
      rec <- rec_vol(nodes[j, ], n, vi, vo)
      vol[j, i] <- rec[1,1]
      lvl[j, i] <- rec[1,2]
    }
  }
  data.frame(
    vol = apply(vol, 1, median),
    lvl = apply(lvl, 1, median))
}

timer <- function(exp, msg = 'time elapsed: ') {
  begin <- Sys.time()
  exp
  end <- Sys.time()
  print(paste0(msg, end - begin))
}


# function to search surface topology for minima

select_rate_pairs <- function(ar, xs, ys) {
  vec <- as.vector(ar)
  ardim <- dim(ar)
  rvec <- rep(1:ardim[1], ardim[2])
  cvec <- 0
  for (i in 1:ardim[2]) {
    cvec <- c(cvec, rep(i, ardim[1]))
  }
  cvec <- cvec[-1]
  cdf <- cumsum(vec)
  xid <- rvec[cdf == min(cdf[cdf > runif(1)])]
  yid <- cvec[cdf == min(cdf[cdf > runif(1)])]
  x <- runif(1, xs[xid], xs[xid+1])
  y <- runif(1, ys[yid], ys[yid+1])
  c(x, y)
}


weight_array <- function(ar, scl = 1) {
  vec <- as.vector(ar)
  wmax <- max(vec)
  p <- (1 - vec) * scl * (1 / length(vec))
  p <- p / sum(p)
  array(p, dim(ar))
}

weight_grid <- function(pos, xy_grid, rec, test, by, scl = 1) {
  val <- 0
  if (test == 'ks') val <- 3
  if (test == 'kp') val <- 5
  if (by == 'lvl') val <- val + 1
  x_dist <- xy_grid[, 1] - pos[1]
  y_dist <- xy_grid[, 2] - pos[2]
  dist <- sqrt(x_dist^2 + y_dist^2)
  dmax <- max(dist) + 1
  wts <- (dmax - dist) / dmax * scl
  sum((wts * rec[, val])) / sum(wts)
}


update_search_array <- function(rec, xs, ys, ar, test, by,
                                d_scl = 1, w_scl = 0.9) {
  xy_grid <- matrix(0, nrow = nrow(rec), ncol = 2)
  for (i in 1:nrow(rec)) {
    x_bin <- 0
    k <- 1
    while (!x_bin) {
      k <- k + 1
      if (rec[i, 1] <= xs[k] & rec[i, 1] > xs[k-1]) x_bin <- k
    }
    y_bin <- 0
    k <- 1
    while (!y_bin) {
      k <- k + 1
      if (rec[i, 2] <= ys[k] & rec[i, 2] > ys[k-1]) y_bin <- k
    }
    xy_grid[i, ] <- c(x_bin, y_bin)
  }
  for (i in 1:nrow(ar)) {
    for (j in 1:ncol(ar)) {
      ar[i, j] <- weight_grid(c(i,j), xy_grid, rec, test, by, d_scl)
    }
  }
  weight_array(ar, w_scl)
}


search_topology <- function(so,
                            x_range,
                            y_range,
                            rec = 0,
                            bt = 10,
                            it = 20,
                            n = 10000,
                            test = 'ks',
                            by = 'lvl',
                            grid = 100,
                            d_scl = 1,
                            w_scl = 3) {
  x_dist <- max(x_range) - min(x_range)
  y_dist <- max(y_range) - min(y_range)
  xs <- seq(
    min(x_range),
    max(x_range),
    x_dist / grid
  )
  ys <- seq(
    min(y_range),
    max(y_range),
    y_dist / grid
  )

  ar <- array(1/(grid^2), c(grid, grid))
  if (length(rec) == 1) vec <- 0
  if (length(rec) > 1) {
    vec <- rec[1, ]
    if (nrow(rec) > 1) {
      for (i in 2:nrow(rec)) {
        vec <- c(vec, rec[i, ])
      }
    }
  }
  begin <- Sys.time()
  for (i in 1:bt) {
    if (length(vec) > 1) {
      ar <- update_search_array(
        matrix(unlist(vec), ncol = 6, byrow = TRUE),
        xs, ys, ar, test, by, d_scl, w_scl)
      present_topology(ar, xs, ys)
    }
    rates <- select_rate_pairs(ar, xs, ys)
    boot <- fit_volumes(so, n, rates[1], rates[2], it)
    vol_gof <- cdf_gof(boot[, 1], so$xsec_area)
    lvl_gof <- cdf_gof(boot[, 2], so$lvl)

    if (length(vec) > 1) vec <- c(vec, rates, vol_gof, lvl_gof)
    if (length(vec) == 1) vec <- c(rates, vol_gof, lvl_gof)
    now <- Sys.time()
    dif <- now - begin
    pace <- dif / i
    rem <- (pace * bt) - dif
    print(paste0('Percent Done: ', round(i/bt*100, 2), '%'))
    print(paste0('Time Elapsed: ', round(as.numeric(dif, units = 'hours'), 2), ' hours'))
    print(paste0('Time Remaining: ', round(as.numeric(rem, units = 'hours'), 2), ' hours'))

  }
  end <- Sys.time()
  print(end - begin)
  matrix(unlist(vec), ncol = 6, byrow = TRUE)
}

record <- search_topology(
  crks, #[crks$creek_name == 'knowles', ],
  c(0,1), c(0,50), rec = rec, bt = 20, by = 'vol', w_scl = 4)
record[record[, 3] < .10, ]
scatter3D(record[,1], record[,2], record[,3],
          ticktype = 'detailed', pch = 20, theta = 240, phi = 40)


rec <- record

best_so_far <- rec
rec <- best_so_far
linked <- rec

save(best_so_far, file = 'bookmark_9-17.rds')
load('bookmark_9-17.rds')

lo_li <- record #9/17, 1000
save(lo_li, file = 'lo_li_9-27_1000.rds')
lpwr_out <- rbind(lpwr_out, record)

record <- best_so_far

best_fit <- record[record[,3] == min(record[,3]), ]

pred <- fit_volumes(crks, 10000, best_fit[1], best_fit[2], 10)
plot(crks$xsec_area %>% emp_cdf, pch = 20, col = get_palette('ocean', .2))
lines(pred$vol %>% emp_cdf, lwd = 3, col = get_palette('charcoal', .6))

plot(crks$lvl %>% emp_cdf)
lines(pred$lvl %>% emp_cdf, lwd = 3, col = get_palette('ocean'))

plot(crks$lvl %>% emp_cdf)
lines(pred1$lvl %>% emp_cdf, lwd = 3, col = get_palette('ocean'))
plot(crks$xsec_area %>% emp_cdf)
lines(pred1$vol %>% emp_cdf, lwd = 3, col = get_palette('ocean'))

linked_bucket <- record
best_so_far <- record
glob_out <- record

png('ri_vs_ro.png', height = 12, width = 12, units = 'cm', res = 300)
plot(crks$ri, crks$ro, pch = 20, col = get_palette('ocean'),
     xlab = 'input rate', ylab = 'output rate')
dev.off()

png('br_kn_xsec_area.png', height = 17, width = 17, units = 'cm', res = 300)
plot(crks$xsec_area %>% emp_cdf, pch = 20, col = get_palette('charcoal'),
     xlab = 'cross-sectional area m2', ylab = 'CDF')
lines(crks$xsec_area %>% emp_cdf, lwd = 2, col = get_palette('ocean'))
dev.off()

png('br_kn_lvl.png', height = 17, width = 17, units = 'cm', res = 300)
plot(crks$lvl %>% emp_cdf, pch = 20, col = get_palette('charcoal'),
     xlab = 'cross-sectional area m2', ylab = 'CDF')
lines(crks$lvl %>% emp_cdf, lwd = 2, col = get_palette('ocean'))
dev.off()


# linked-bucket model
# pass removals downstream
# to account for fluvial fines and gravels

crks$rdp <- crks$dp / max(crks$dp)

record_volume <- function(node, n, vi, vo,
                          arrivals = 0) {
  departures <- 0
  if (length(arrivals) != 1) {
    rolls <- runif(nrow(arrivals))
    if (sum(rolls > node$rdp) >= 1) {
      departures <- arrivals[rolls > node$rdp, ]
      arrivals <- arrivals[rolls <= node$rdp, ]
    }
    arrivals$lvl <- arrivals$vol / as.numeric(node$valley_width)
  }
  t <- 0
  age <- 0
  vol <- 0
  k <- 0
  while (t < n) {
    k <- k + 1
    t <- t + rexp(1, node$ri)
    age[k] <- t
    vol[k] <- rexp(1, vi)
  }
  df <- data.frame(
    t = age,
    r = 'input',
    vol = vol,
    lvl = vol / as.numeric(node$valley_width)
  )
  if (length(arrivals) != 1) df <- rbind(df, arrivals)
  t <- 0
  age <- 0
  vol <- 0
  k <- 0
  while (t < n) {
    k <- k + 1
    t <- t + rexp(1, node$ro)
    age[k] <- t
    vol[k] <- rexp(1, vo)
  }
  dep <- data.frame(
    t = age,
    r = 'output',
    vol = vol,
    lvl = vol / as.numeric(node$valley_width))
  if (length(departures) == 1) departures <- dep
  if (length(departures) != 1) departures <- rbind(departures, dep)
  departures$r <- 'input'
  df <- rbind(df, dep)
  df <- df[order(df$t),]
  vol <- 0
  lvl <- 0
  for (i in 1:nrow(df)) {
    if (df$r[i] == 'input') {
      vol <- vol + df$vol[i]
      lvl <- lvl + df$lvl[i]
    }
    if (df$r[i] == 'output') {
      if (df$vol[i] >= vol) {
        vol <- 0
      }
      if (df$vol[i] < vol) {
        vol <- vol - df$vol[i]
      }
      if (df$lvl[i] >= lvl) {
        lvl <- 0
      }
      if (df$lvl[i] < lvl) {
        lvl <- lvl - df$lvl[i]
      }
    }
  }
  res <- data.frame(vol = vol, lvl = lvl)
  list(res, departures)
}

# function to fit volume rate to observed volumes

fit_volumes1 <- function(nodes, n, vi, vo, it = 10) {
  vol <- array(0, c(nrow(nodes), it))
  lvl <- array(0, c(nrow(nodes), it))
  nodes <- nodes[order(nodes$ToMouth_km, decreasing = TRUE), ]
  arrivals <- 0
  for (i in 1:it) {
    for (j in 1:nrow(nodes)) {
      if (length(arrivals) == 1) res <- record_volume(nodes[j, ], n, vi, vo)
      if (length(arrivals) > 1) {
        res <- record_volume(nodes[j, ], n, vi, vo, arrivals)
      }
      rec <- res[[1]]
      arrivals <- res[[2]]
      vol[j, i] <- rec[1,1]
      lvl[j, i] <- rec[1,2]
    }
  }
  data.frame(
    vol = apply(vol, 1, median),
    lvl = apply(lvl, 1, median))
}


# function to visualize topology

present_topology <- function(ar, xs, ys) {
  xvec <- rep(1:nrow(ar), ncol(ar))
  yvec <- c(
    matrix(rep(1:ncol(ar), nrow(ar)), ncol = ncol(ar), byrow = TRUE)
    )
  vec <- 1 - (c(ar) / max(ar))
  print(scatter3D(xs[xvec], ys[yvec], vec, pch = 20, ticktype = 'detailed',
                  phi = 40, theta = 320))
  return()
}







bfs_vol_lwr <- data.frame(ri = csq[,2],
                         ro = csq[,1],
                         ks = gofs[,1],
                         kp = gofs[,2])
save(bfs_vol_lwr, file = 'bfs_vol_lwr.rds')

#png('bfs_vol_lwr.png', height = 30, width = 30, units = 'cm', res = 300)
scatter3D(bfs_vol_lwr[,1], bfs_vol_lwr[,2], bfs_vol_lwr[,3], ticktype = 'detailed', pch = 20,
          xlab = 'output rate', ylab = 'input rate', zlab = 'K-S stat',
          phi = 35, theta = 150)
#dev.off()


dim(crks)
scatter3D(bfs_io[,1], bfs_io[,2], bfs_io[,3], pch = 20)
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

ksmin <- bfs_ia[bfs_ia$ks == min(bfs_ia$ks),]
kpmin <- bfs_ia[bfs_ia$kp == min(bfs_ia$kp),]

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


sz <- length(dfmns)
bts <- 100
plot(emp_cdf(dfmns),
     xlab = 'transit time', ylab = 'cdf')
for (i in 1:nrow(ksmin)) {
  lines(emp_cdf(accumulate(sz, ksmin[i,1], ksmin[i,2], ia_ar, bts)),
        lwd = 2, col = get_palette('crimson', .33))
}
for (i in 1:nrow(kpmin)) {
  lines(emp_cdf(accumulate(sz, kpmin[i,1], kpmin[i,2], ia_ar, bts)),
        lwd = 2, col = get_palette('ocean', .33))
}
legend('bottomright', legend = c('observed', 'K-S fit', 'Kuiper fit'),
       fill = get_palette(c('charcoal', 'crimson', 'ocean'), .7))










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

dfmns <- read.csv('dfmns.csv')
dfmns <- dfmns$age %>% unlist

library(plot3D)
rec <- fread('/home/crumplecup/projects/reservoir/df_boot.csv')
med <- apply(rec, 2, sort)
med <- apply(med, 1, median)
png('test.png', height = 15, width = 18, units ='cm', res = 300)
plot(emp_cdf(med),
     pch = 20, col = get_palette('ocean'))
lines(emp_cdf(dfmns), lwd = 3, col = get_palette('forest'))
points(emp_cdf(dfmns), pch = 20, col = get_palette('forest'))
lines(emp_cdf(dfmn), lwd = 3, col = get_palette('gold'))
points(emp_cdf(dfmn), pch = 20, col = get_palette('gold'))
dev.off()

rec <- fread('/home/crumplecup/projects/reservoir/dfio_fit2.csv')
nrow(rec)

scatter3D(rec$input, rec$output, rec$ks, pch = 20, phi = 30, theta = 50,
          xlab = "input", ylab = 'output', zlab = 'fit',
          ticktype = 'detailed')
rec[rec$ks < min(rec$ks)*1.1, ]


sub2 <- rec[rec$ff_ks > 0 & rec$ff_n < 200, ]
scatter3D(sub2$fup, sub2$ff_n, sub2$ff_ks, pch = 20, phi = 30, theta = 50,
          xlab = "uptake", ylab = 'samples', zlab = 'fit',
          ticktype = 'detailed')
sub2 <- rec[rec$fg_ks > 0 & rec$fg_n < 200, ]
scatter3D(sub2$gup, sub2$fg_n, sub2$fg_ks, pch = 20, phi = 10, theta = 30,
          xlab = "uptake", ylab = 'samples', zlab = 'fit',
          ticktype = 'detailed')

sub2 <- rec[rec$fg_ks > 0 & rec$fg_n < 90.5 & rec$fg_n > 89.5, ]
scatter3D(sub2$gup, sub2$fg_n, sub2$fg_ks, pch = 20, phi = 10, theta = 30,
          xlab = "uptake", ylab = 'samples', zlab = 'fit',
          ticktype = 'detailed')
sub2[sub2$gup == max(sub2$gup), ]

sub2 <- rec[rec$ff_ks > 0 & rec$ff_n < 56.5 & rec$ff_n > 55.5, ]
scatter3D(sub2$fup, sub2$ff_n, sub2$ff_ks, pch = 20, phi = 10, theta = 30,
          xlab = "uptake", ylab = 'samples', zlab = 'fit',
          ticktype = 'detailed')
sub2[sub2$fup == max(sub2$fup), ]

sub <- sub2[sub2$ff_ks < .10, ]
scatter3D(sub$up, sub$fro, sub$ff_n, pch = 20, phi = 30, theta = 200,
          xlab = "uptake", ylab = 'ff out', zlab = 'k-s value',
          ticktype = 'detailed')
rec[rec$fg_ks < .1 & rec$fg_ks > 0, ]

plot(ffsub$fro, ffsub$ff_ks, pch = 20, col = get_palette('gold'))
points(sub2$fro, sub2$ff_ks, pch = 20, col = get_palette('ocean'))


# sub <- rec[rec$df_ks < .132, ]
# scatter3D(sub$ri, sub$ro, sub$df_ks, pch = 20, phi = 30, theta = 230,
#           xlab = "rate in", ylab = 'rate out', zlab = 'k-s value',
#           ticktype = "detailed")

# sub <- rec[rec$ff_ks < .5, ]
# scatter3D(sub$ri, sub$ro, sub$ff_ks, pch = 20, phi = 30, theta = 230,
#           xlab = "rate in", ylab = 'rate out', zlab = 'k-s value',
#           ticktype = "detailed")


