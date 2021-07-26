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


gravel_fit_n <- function(n, capture, storage, turnover) {
  mat <- matrix(0, n, 4)
  mat <- apply(mat, 1, function(x) gof(gravel_synth(capture, storage, turnover), gravels))
  apply(mat, 1, mean)
}

fines_fit_n <- function(n, capture, storage, turnover) {
  mat <- matrix(0, n, 4)
  mat <- apply(mat, 1, function(x) gof(fines_synth(capture, storage, turnover), fines))
  apply(mat, 1, mean)
}


gravel_fit_n(10, .12, .12, 208)
fines_fit_n(10, .12, .12, 208)

gravel_fit <- function(batch = 10, n = 10, turnover, max_cap = 1, max_stor = 1, fines = F) {
  capture_rates <- runif(batch, 0, max_cap)
  storage_rates <- runif(batch, 0, max_stor)
  turnovers <- runif(batch, 175, 350)
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
    ad = mat[ , 1],
    ch = mat[ , 2],
    kp = mat[ , 3],
    ks = mat[ , 4])
}

# function to split observations of two predictors into bins
# returns the mean sd and count of fit for each bin
bin_dual_stat <- function(pred1, pred2, fit, bins = 10) {
  # vector of means, sds and counts
  mns <- 0
  sds <- 0
  ns <- 0
  lwr <- 0
  upr <- 0
  # divide range of pred into bins number of steps
  rng1 <- seq(min(pred1), max(pred1), (max(pred1) - min(pred1)) / bins)
  rng2 <- seq(min(pred2), max(pred2), (max(pred2) - min(pred2)) / bins)
  print(length(rng1))
  print(length(rng2))
  for (i in 1:bins) {
    for (j in 1:bins) {

    }
    # select fits where pred in is range
    if (i == 1 & j == 1) {
      bin <- fit[pred1 <= rng1[i] & pred2 <= rng2[j]]
    } else if (i == 1 & j > 1) {
      bin <- fit[pred1 <= rng1[i] &
                   pred2 > rng2[j-1] & pred2 <= rng2[j]]
    } else if (i > 1 & j == 1) {
      bin <- fit[pred1 > rng1[i-1] & pred1 <= rng1[i] &
                   pred2 <= rng2[j]]
    } else {
      bin <- fit[pred1 > rng1[i-1] & pred1 <= rng1[i] &
                   pred2 > rng2[j-1] & pred2 <= rng2[j]]
    }
    mns[i] <- mean(bin)
    sds[i] <- sd(bin)
    ns[i] <- length(bin)
    lwr[i] <- mns[i] - 1.96 * (sds[i] / sqrt(ns[i]))
    upr[i] <- mns[i] + 1.96 * (sds[i] / sqrt(ns[i]))
  }
  data.frame(mns = mns, sds = sds, ns = ns, upr = upr, lwr = lwr, rng1 = rng1[1:length(mns)], rng2 = rng2[1:length(mns)])
}

rec[rec$ch == min(rec$ch), ]
gr_stats <- bin_dual_stat(rec$capture, rec$storage, rec$ks, 20)[-1, ]
gr_stats$rng1[gr_stats$lwr <= gr_stats$upr[gr_stats$mns == min(gr_stats$mns)]]
gr_stats$rng2[gr_stats$lwr <= gr_stats$upr[gr_stats$mns == min(gr_stats$mns)]]
plot3D::points3D(gr_stats$rng1, gr_stats$rng2, gr_stats$mns)

gravel_fit(10, 10, 208, 0.2, 0.2)

rec <- gravel_fit(10, 200, 318, 0.5, 1)
rec <- gravel_fit(10, 200, 318, 1, 1, fines = T)
dur <- as.difftime(3, units = 'hours')
begin <- Sys.time()
end <- begin + dur
while (Sys.time() < end) {
  rec <- rbind(rec, gravel_fit(10, 200, 318, 1, 1, fines = T))
  save(rec, file = 'fn_200.rds')
}

rec[rec$ks == min(rec$ks), ]

plot3D::points3D(rec$capture, rec$storage, rec$ad, ticktype = 'detailed', pch = 20,
                 phi = 20, theta = 210,
                 xlab = 'capture rate', ylab = 'storage rate')

plot3D::points3D(rec$capture, rec$storage, rec$ch, ticktype = 'detailed', pch = 20,
                 phi = 20, theta = 340,
                 xlab = 'capture rate', ylab = 'storage rate')

plot3D::points3D(rec$capture, rec$storage, rec$kp, ticktype = 'detailed', pch = 20,
                 phi = 20, theta = 310,
                 xlab = 'capture rate', ylab = 'storage rate')

plot3D::points3D(rec$capture, rec$storage, rec$ks, ticktype = 'detailed', pch = 20,
                 phi = 20, theta = 310,
                 xlab = 'capture rate', ylab = 'storage rate')

magicaxis::magplot(rec$storage, rec$ks, pch = 20, col = get_palette('ocean'))

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










