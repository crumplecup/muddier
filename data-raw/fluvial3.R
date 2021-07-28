library(magrittr)
library(parallel)
options(mc.cores=parallel::detectCores())
setwd('/media/erik/catacomb/research')

debris_flows <- charcoal$mn[charcoal$facies == 'DF'] + 50
gravels <- charcoal$mn[charcoal$facies == 'FG'] + 50
fines <- charcoal$mn[charcoal$facies == 'FF'] + 50

# debris-flow pmf and index
index <- char_pmfs %>% rownames %>% as.numeric %>% rev
dfc <- emp_cdf(charcoal$mn[charcoal$facies == 'DF'])
dfi <- dfc[ , 1]
dfi <- dfi + 50
dfp <- to_pmf(dfc[ , 2])

# gravels pmf and index
grc <- emp_cdf(gravels)
gri <- grc[ , 1]
grp <- to_pmf(grc[ , 2])


fines_synth <- function(capture, storage, turnover) {
  ages <- sample(
    c(sample(dfi, length(dfi), replace = TRUE, prob = dfp),
      sample(gri, length(gri), replace = TRUE, prob = grp))
    )
  storage_pmf <- fish(storage, index / turnover, 1)
  storage_pmf <- storage_pmf / sum(storage_pmf)
  for (i in 1:length(dfi)) {
    if (runif(1) <= capture) {
      ages[i] <- ages[i] + sample(index, 1, prob = storage_pmf)
    }
  }
  ages
}

gravel_synth <- function(capture, storage, turnover) {
  ages <- sample(dfi, length(dfi), replace = TRUE, prob = dfp)
  storage_pmf <- fish(storage, index / turnover, 1)
  storage_pmf <- storage_pmf / sum(storage_pmf)
  for (i in 1:length(dfi)) {
    if (runif(1) <= capture) {
      ages[i] <- ages[i] + sample(index, 1, prob = storage_pmf)
    }
  }
  ages
}

# poisson distribution
# given a rate, time and number of events
# returns p(k) distribution
fish <- function(rt, t, k) {
  (rt*t)^k*exp(-rt*t) / factorial(k)
}

plot(emp_cdf(gravel_synth(.12,.12,318)))


gof <- function(synth, obs) {
  vals <- sort(unique(c(synth, obs)))
  kobs <- sort(c(synth, obs))
  lnx <- length(synth)
  lny <- length(obs)
  k <-  lnx + lny
  c1 <- 0
  c1l <- 0
  c2 <- 0
  c2l <- 0
  k1 <- 0
  chi <- 0
  for (i in seq_along(vals)) {
    c1l[i] <- length(synth[synth <= vals[i]])
    c1[i] <-  c1l[i] / length(synth)
    c2l[i] <- length(obs[obs <= vals[i]])
    c2[i] <-  c2l[i] / length(obs)
    chi[i] <- (c1l[i] - c2l[i])^2 / c2l[i]
    k1[i] <- length(kobs[kobs <= vals[i]]) / length(kobs)
  }
  ks <- max(abs(c1 - c2))
  kp1 <- max(c1 - c2)
  kp2 <- max(c2 - c1)
  kp <- kp1 + kp2
  ch <- sum(chi[chi != Inf])
  adi <- 0
  for (i in 1:(k-1)) {
    xl <- length(synth[synth <= kobs[i]])
    adi[i] <- (lnx * xl - lnx * i)^2 / (i * (k - i))
  }
  ad <- (1 / (lnx * lny)) * sum(adi)
  return(c(ad, ch, kp, ks))
}

test_fit <- function(fits) {
  hits <- c(0,0,0,0, 0,0,0,0)
  # anderson-darling sig thresholds
  if (fits[1] < 340) hits[1] <- 1
  if (fits[1] < 250) hits[5] <- 1
  # chi-squared sig thresholds
  if (fits[2] < 300) hits[2] <- 1
  if (fits[2] < 100) hits[6] <- 1
  # kuiper sig thresholds
  if (fits[3] < 0.12) hits[3] <- 1
  if (fits[3] < 0.09) hits[7] <- 1
  # kolmogorov-smirnov sig thresholds
  if (fits[4] < 0.10) hits[4] <- 1
  if (fits[4] < 0.05) hits[8] <- 1
  hits
}

gravel_fit_n <- function(n, capture, storage, turnover) {
  mat <- matrix(0, n, 4)
  mat <- apply(mat, 1, function(x) gof(gravel_synth(capture, storage, turnover), gravels))
  res <- matrix(0, nrow = n, ncol = 8)
  for (i in 1:n) {
    res[i, ] <- test_fit(mat[, i])
  }
  apply(res, 2, function(x) sum(x) / length(x))
}

gravel_fit_n(10, .12, .12, 208)

fines_fit_n <- function(n, capture, storage, turnover) {
  mat <- matrix(0, n, 4)
  mat <- apply(mat, 1, function(x) gof(fines_synth(capture, storage, turnover), fines))
  res <- matrix(0, nrow = n, ncol = 8)
  for (i in 1:n) {
    res[i, ] <- test_fit(mat[, i])
  }
  apply(res, 2, function(x) sum(x) / length(x))
}

fines_fit_n(10, .12, .12, 208)

gravel_fit <- function(batch = 10, n = 10, max_cap = 1, max_stor = 1, fines = F) {
  capture_rates <- runif(batch, 0, max_cap)
  storage_rates <- runif(batch, 0, max_stor)
  turnovers <- runif(batch, 50, 1500)
  mat <- 0
  if (fines) {
    mat <- t(parallel::mcmapply(function(a,b,c,d,e) fines_fit_n(a, b, c, d),
                                b = capture_rates, c = storage_rates, d = turnovers,
                                MoreArgs = list(a = n)))
  } else {
    mat <- t(parallel::mcmapply(function(a,b,c,d,e) gravel_fit_n(a, b, c, d),
                                b = capture_rates, c = storage_rates, d = turnovers,
                                MoreArgs = list(a = n)))
  }
  data.frame(
    capture = capture_rates,
    storage = storage_rates,
    turnover = turnovers,
    ad_a = mat[ , 1],
    ch_a = mat[ , 2],
    kp_a = mat[ , 3],
    ks_a = mat[ , 4],
    ad_b = mat[ , 5],
    ch_b = mat[ , 6],
    kp_b = mat[ , 7],
    ks_b = mat[ , 8])
}

# function to split observations of two predictors into bins
# returns the mean sd and count of fit for each bin
bin_dual_stat <- function(pred1, pred2, fit, bins = 10) {
  # vector of means, sds and counts
  mns <- 0
  sds <- 0
  ns <- 0
  idx1 <- 0
  idx2 <- 0

  # divide range of pred into bins number of steps
  rng1 <- seq(min(pred1), max(pred1), (max(pred1) - min(pred1)) / bins)
  rng2 <- seq(min(pred2), max(pred2), (max(pred2) - min(pred2)) / bins)
  for (i in 1:bins) {
    if (i == 1) {
      bin <- fit[pred1 <= rng1[i]]
    } else {
      bin <- fit[pred1 > rng1[i-1] & pred1 <= rng1 [i]]
    }
    for (j in 1:bins) {
      if (j == 1) {
        jbin <- bin[pred2 <= rng2[j]]
      } else {
        jbin <- bin[pred2 > rng2[j-1] & pred2 <= rng2[j]]
      }
      mns <- c(mns, mean(jbin[!is.na(jbin)]))
      sds <- c(sds, sd(jbin[!is.na(jbin)]))
      ns <- c(ns, length(jbin[!is.na(jbin)]))
      idx1 <- c(idx1, rng1[i])
      idx2 <- c(idx2, rng2[j])
    }
  }
  mns <- mns[-1]
  sds <- sds[-1]
  ns <- ns[-1]
  idx1 <- idx1[-1]
  idx2 <- idx2[-1]

  lwr <- mns - 1.96 * (sds / sqrt(ns))
  upr <- mns + 1.96 * (sds / sqrt(ns))

  data.frame(mns = mns, sds = sds, ns = ns, upr = upr, lwr = lwr, rng1 = idx1, rng2 = idx2)
}


gravel_fit(10, 10, 1, 1)


load('/home/erik/output/gr_1000.rds')
rc <- rec
load('/home/erik/output/gr_1000a.rds')
rc <- rbind(rc, rec)
load('gr_1000.rds')
rc <- rbind(rc, rec)
rec[rec$ch == min(rec$ch), ]
(gr_stats <- bin_dual_stat(rec$capture, rec$storage, rec$ks, 20))
gr_min <- min(gr_stats$mns[!is.na(gr_stats$mns) & gr_stats$mns > 0])
range(gr_stats$rng1[gr_stats$lwr <= gr_stats$upr[gr_stats$mns == gr_min & !is.na(gr_stats$upr)]
              & !is.na(gr_stats$lwr) & gr_stats$lwr > 0])
range(gr_stats$rng2[gr_stats$lwr <= gr_stats$upr[gr_stats$mns == gr_min & !is.na(gr_stats$upr)]
              & !is.na(gr_stats$lwr) & gr_stats$lwr > 0])
plot3D::points3D(gr_stats$rng1, gr_stats$rng2, gr_stats$mns)

load('/home/erik/output/fn_1000.rds')
rc <- rec
load('/home/erik/output/fn_200.rds')
rc <- rbind(rc, rec)
load('fn_200.rds')
rc <- rbind(rc, rec)
rec[rec$ks == min(rec$ks), ]
(gr_stats <- bin_dual_stat(rec$capture, rec$storage, rec$ks, 20)[-1, ])
gr_stats$rng1[gr_stats$lwr <= gr_stats$upr[gr_stats$mns == min(gr_stats$mns)]]
gr_stats$rng2[gr_stats$lwr <= gr_stats$upr[gr_stats$mns == min(gr_stats$mns)]]
plot3D::points3D(gr_stats$rng1, gr_stats$rng2, gr_stats$mns)



rec <- gravel_fit(10, 200, 1, 1)
# rec <- gravel_fit(10, 200, 318, 1, 1, fines = T)
# load('fn_200.rds')
dur <- as.difftime(3, units = 'hours')
begin <- Sys.time()
end <- begin + dur
while (Sys.time() < end) {
  rec <- rbind(rec, gravel_fit(10, 200, 1, 1))
  save(rec, file = 'gr_200_hits.rds')
}

rec[rec$ks == min(rec$ks), ]
rec[rec$kp == min(rec$kp), ]

rec[rec$ad_a == max(rec$ad_a), ]
plot3D::points3D(rec$capture, rec$storage, rec$ad_a, ticktype = 'detailed', pch = 20,
                 phi = 40, theta = 30,
                 xlab = 'capture rate', ylab = 'storage rate', zlab = 'ad alpha')

rec[rec$ch_a == max(rec$ch_a), ]
plot3D::points3D(rec$capture, rec$storage, rec$ch_a, ticktype = 'detailed', pch = 20,
                 phi = 20, theta = 290,
                 xlab = 'capture rate', ylab = 'storage rate')

rec[rec$kp_a == max(rec$kp_a), ]
plot3D::points3D(rec$capture, rec$storage, rec$kp_a, ticktype = 'detailed', pch = 20,
                 phi = 30, theta = 310,
                 xlab = 'capture rate', ylab = 'storage rate')

rec[rec$ks_a == max(rec$ks_a), ]
plot3D::points3D(rec$capture, rec$storage, rec$ks_a, ticktype = 'detailed', pch = 20,
                 phi = 40, theta = 310,
                 xlab = 'capture rate', ylab = 'storage rate')

magicaxis::magplot(rec$turnover, rec$kp_b, pch = 20, col = get_palette('ocean'),
                   xlab = 'turnover period [y]', ylab = 'k-s fit')
abline(v = 191, col = get_palette('violet', .9), lwd = 2)
abline(v = 208, col = get_palette('crimson', .9), lwd = 2)
abline(v = 293, col = get_palette('gold', .9), lwd = 2)
abline(v = 318, col = get_palette('ocean', .9), lwd = 2)

# load('gr_1000.rds')
sub <- rec
subad <- sub[sub$turnover > 300 , ]
subch <- sub[sub$turnover > 290 & sub$turnover < 300 , ]
subkp <- sub[sub$turnover < 200 , ]
subks <- sub[sub$turnover > 200 & sub$turnover < 215 , ]
magicaxis::magplot(sub$turnover, sub$ks, pch = 20, col = get_palette('ocean'))

plot3D::points3D(subad$capture, subad$storage, subad$ks, ticktype = 'detailed', pch = 20,
                 phi = 20, theta = 310,
                 xlab = 'capture rate', ylab = 'storage rate')
plot3D::points3D(subch$capture, subch$storage, subch$ks, ticktype = 'detailed', pch = 20,
                 phi = 20, theta = 310,
                 xlab = 'capture rate', ylab = 'storage rate')
plot3D::points3D(subkp$capture, subkp$storage, subkp$ks, ticktype = 'detailed', pch = 20,
                 phi = 20, theta = 310,
                 xlab = 'capture rate', ylab = 'storage rate')
plot3D::points3D(subks$capture, subks$storage, subks$ks, ticktype = 'detailed', pch = 20,
                 phi = 20, theta = 310,
                 xlab = 'capture rate', ylab = 'storage rate')

magicaxis::magplot(emp_cdf(charcoal$mn[charcoal$facies == 'FF'] + 50), pch = 20, col = get_palette('coral'),
                   xlab = 'charcoal age', ylab = 'CDF of samples', xlim = c(0, 17000))
points(emp_cdf(debris_flows), pch = 20, col = get_palette('hardwood'))
points(emp_cdf(charcoal$mn[charcoal$facies == 'FG'] + 50), pch = 20, col = get_palette('charcoal'))
lines(emp_cdf(gravel_synth(.07, .04, 191)))


n <- 10
gravels <- data.frame(
  capture = runif(n, 0, 0.5),
  storage = runif(n, 0, 0.5)
)

res <- list()
for (i in 1:nrow(gravels)) {
  res[[i]] <- gravel_synth(gravels[i, 1], gravels[i, 2], 208)
}
lapply(gravels, function(x) gravel_synth(x[1], x[2], 208))


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
fines_synth <- function(storage = ffep, capture = fep, n = 10) {
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


tp_ad <- 318
tp_ch <- 293
tp_kp <- 191
tp_ks <- 208

t <- seq(0, 37, length.out = 4430)


pex_ad <- (fish(pf_ad^trap, t, 1)/sum(fish(pf_ad^trap,t,1)))
pex_ch <- (fish(pf_ch^trap, t, 1)/sum(fish(pf_ch^trap,t,1)))
pex_kp <- (fish(pf_kp^trap, t, 1)/sum(fish(pf_kp^trap,t,1)))
pex_ks <- (fish(pf_ks^trap, t, 1)/sum(fish(pf_ks^trap,t,1)))



fns_ks <- 0
fns_kp <- 0
for (i in 1:10) {
  fns_kp <- c(fns_kp, fines_synth(.301, .04, 191))
  fns_ks <- c(fns_ks, fines_synth(.222, .054, 208))
}
fns_kp <- fns_kp[-1]
fns_ks <- fns_ks[-1]
lines(emp_cdf(fns_kp))


df_trans_ks <- data.table::fread('/home/erik/output/transits_cdf_gr_kp.csv')
gr_trans_ks <- data.table::fread('/home/erik/output/transits_cdf_gr_ks.csv')
gr_trans_ks <- data.table::fread('/home/erik/output/gravels_transits_ks.csv')
gr_trans_ks <- data.table::fread('/home/erik/output/gravels_transits_kp.csv')
gr_trans_kp <- data.table::fread('/home/erik/output/transits_cdf_gr_kp.csv')
gr_trans_ch <- data.table::fread('/home/erik/output/transits_cdf_gr_ch.csv')
gr_trans_ks <- data.table::fread('/home/erik/output/transits_cdf_fn_ks.csv')
gr_trans_kp <- data.table::fread('/home/erik/output/transits_cdf_fn_kp.csv')
df_trans_ks <- data.table::fread('/home/erik/output/debris_flow_deposits_cdf.csv')


png('fines_fit.png', height = 17, width = 21, units = 'cm', res = 300)
magicaxis::magplot(emp_cdf(charcoal$mn[charcoal$facies == 'FF'] + 50), pch = 20, col = get_palette('coral'),
                   xlab = 'Charcoal Age [y]', ylab = 'CDF of samples', xlim = c(0, 17000))
points(emp_cdf(debris_flows), pch = 20, col = get_palette('hardwood'))
points(emp_cdf(charcoal$mn[charcoal$facies == 'FG'] + 50), pch = 20, col = get_palette('charcoal'))

lines(emp_cdf(fns_ks), lwd = 2.5, col = get_palette('crimson', .7))
lines(emp_cdf(fns_kp), lwd = 2.5, col = get_palette('violet', .7))
# lines(0:20000, cumsum(unlist(gr_trans_ch)), lwd = 2.5, col = get_palette('gold', .7))
legend('bottomright', legend = c('Debris-Flows', 'Gravels', 'Fines', 'Kolmogorov-Smirnov', 'Kuiper'),
       lty = c(NA, NA, NA, 1, 1), lwd = c(NA, NA, NA, 2, 2),
       pch = c(20, 20, 20, NA, NA), col = get_palette(c('hardwood', 'charcoal', 'coral', 'crimson', 'violet'), .9))
dev.off()







