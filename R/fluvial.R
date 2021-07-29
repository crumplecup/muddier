#' fines_synth
#'
#' simulate fines ages
#'
#' @param capture the capture rate 0-1
#' @param storage the storage rate 0-1
#' @param turnover the turnover period
#' @param storage_index is the age value index for storage
#' @param source_index is the age value index corresponding to the PMF source_prob
#' @param source_prob is the PMF of debris-flow deposit ages
#' @param gravel_index is the age value index corresponding to the PMF gravel_prob
#' @param gravel_prob is the PMF of gravel ages
#' @return vector of synthetic fines ages
#' @export
fines_synth <- function(capture, storage, turnover,
                        storage_index, source_index, source_prob,
                        gravel_index, gravel_prob) {
  ages <- sample(
    c(sample(source_index, length(source_index), replace = TRUE, prob = source_prob),
      sample(gravel_index, length(gravel_index), replace = TRUE, prob = gravel_prob))
  )

  storage_pmf <- fish(storage, storage_index / turnover, 1)
  storage_pmf <- storage_pmf / sum(storage_pmf)
  for (i in 1:length(source_index)) {
    if (runif(1) <= capture) {
      ages[i] <- ages[i] + sample(storage_index, 1, prob = storage_pmf)
    }
  }
  ages
}

#' gravel_synth
#'
#' simulate gravels ages
#'
#' @param capture the capture rate 0-1
#' @param storage the storage rate 0-1
#' @param turnover the turnover period
#' @param storage_index is the age value index for storage
#' @param source_index is the age value index corresponding to the PMF source_prob
#' @param source_prob is the PMF of debris-flow deposit ages
#' @return vector of synthetic fines ages
#' @export
gravel_synth <- function(capture, storage, turnover,
                         storage_index, source_index, source_prob) {
  ages <- sample(source_index, length(source_index), replace = TRUE, prob = source_prob)
  storage_pmf <- fish(storage, storage_index / turnover, 1)
  storage_pmf <- storage_pmf / sum(storage_pmf)
  for (i in 1:length(source_index)) {
    if (runif(1) <= capture) {
      ages[i] <- ages[i] + sample(storage_index, 1, prob = storage_pmf)
    }
  }
  ages
}

#' fish
#'
#' Given a rate, time and number of events,
#' returns a p(k) Poisson distribution.
#'
#' @param rt is the rate.
#' @param t is the time interval.
#' @param k is the number of events.
#' @return a p(k) Poisson distribution, the probability of k events in time t.
#' @export
fish <- function(rt, t, k) {
  (rt*t)^k*exp(-rt*t) / factorial(k)
}


#' gof
#'
#' Goodness-of-fit statistics for fluvial deposits.
#'
#' @param synth Vector of snythetic deposit ages.
#' @param obs Vector of observed deposit ages.
#' @return Vector of fits = c(Anderson-Darling, Chi-Squared, Kuiper, Kolmogorov-Smirnov)
#' @export
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

#' test_fit
#'
#' Hit statistics on how often a vector of fits passes a set of thresholds.
#'
#' @param fits Vector of fits from gof().
#' @param ad_a Less stringent threshold for Anderson-Darling test.
#' @param ad_b More stringent threshold for Anderson-Darling test.
#' @param ch_a Less stringent threshold for Chi-squared test.
#' @param ch_b More stringent threshold for Chi-squared test.
#' @param kp_a Less stringent threshold for Kuiper test.
#' @param kp_b More stringent threshold for Kuiper test.
#' @param ks_a Less stringent threshold for Kolmogorov-Smirnov test.
#' @param ks_b More stringent threshold for Kolmogorov-Smirnov test.
#' @return Vector of hit % = c(ad_a, ch_a, kp_a, ks_a, ad_b, ch_b, kp_b, ks_b)
#' @export
test_fit <- function(fits,
                     # ad_a = 340,
                     # ad_b = 280,
                     ad_a = 240,
                     ad_b = 180,
                     # ch_a = 200,
                     # ch_b = 150,
                     ch_a = 1950,
                     ch_b = 1700,
                     kp_a = 0.14,
                     kp_b = 0.10,
                     ks_a = 0.09,
                     ks_b = 0.05) {
  hits <- c(0,0,0,0, 0,0,0,0)
  # anderson-darling sig thresholds
  if (fits[1] < ad_a) hits[1] <- 1
  if (fits[1] < ad_b) hits[5] <- 1
  # chi-squared sig thresholds
  if (fits[2] < ch_a) hits[2] <- 1
  if (fits[2] < ch_b) hits[6] <- 1
  # kuiper sig thresholds
  if (fits[3] < kp_a) hits[3] <- 1
  if (fits[3] < kp_b) hits[7] <- 1
  # kolmogorov-smirnov sig thresholds
  if (fits[4] < ks_a) hits[4] <- 1
  if (fits[4] < ks_b) hits[8] <- 1
  hits
}

#' gravel_fit_n
#'
#' Fit multiple runs of gravel_synth() to observed deposits.
#'
#' @param n The number of runs to simulate.
#' @param capture The capture rate, ranging 0-1.
#' @param storage The storage rate, ranging 0-1.
#' @param turnover The turnover period, in years.
#' @param storage_index is the age value index for storage
#' @param source_index is the age value index corresponding to the PMF source_prob
#' @param source_prob is the PMF of debris-flow deposit ages
#' @param gravels A vector of observed gravel ages.
#' @return mean hit % of test thresholds, using test_fit()
#' @export
#' @seealso gravel_synth
#' @seealso test_fit
gravel_fit_n <- function(n, capture, storage, turnover, storage_index,
                         source_index, source_prob, gravels) {
  mat <- matrix(0, n, 4)
  mat <- apply(mat, 1, function(x) gof(gravel_synth(capture, storage,
                                                    turnover, storage_index,
                                                    source_index, source_prob),
                                       gravels))
  res <- matrix(0, nrow = n, ncol = 8)
  for (i in 1:n) {
    res[i, ] <- test_fit(mat[, i])
  }
  apply(res, 2, function(x) sum(x) / length(x))
}



#' fines_fit_n
#'
#' Fit multiple runs of fines_synth() to observed deposits.
#'
#' @param n The number of runs to simulate.
#' @param capture The capture rate, ranging 0-1.
#' @param storage The storage rate, ranging 0-1.
#' @param turnover The turnover period, in years.
#' @param storage_index is the age value index for storage
#' @param source_index is the age value index corresponding to the PMF source_prob
#' @param source_prob is the PMF of debris-flow deposit ages
#' @param gravel_index is the age value index corresponding to the PMF gravel_prob
#' @param gravel_prob is the PMF of gravel ages
#' @param fines A vector of observed gravel ages.
#' @return mean hit % of test thresholds, using test_fit()
#' @export
#' @seealso fines_synth
#' @seealso test_fit
fines_fit_n <- function(n, capture, storage, turnover, storage_index,
                        source_index, source_prob,
                        gravel_index, gravel_prob, fines) {
  mat <- matrix(0, n, 4)
  mat <- apply(mat, 1, function(x) gof(fines_synth(capture, storage, turnover, storage_index,
                                                   source_index, source_prob,
                                                   gravel_index, gravel_prob),
                                       fines))
  res <- matrix(0, nrow = n, ncol = 8)
  for (i in 1:n) {
    res[i, ] <- test_fit(mat[, i])
  }
  apply(res, 2, function(x) sum(x) / length(x))
}


#' fluvial_fit
#'
#' Fit multiple fluvial deposit runs to observed deposits.
#'
#' @param batch The number of times to randomly select model parameters.
#' @param n The number of runs to simulate at selected paramaters.
#' @param min_cap = The minimum capture rate to consider, ranging 0-1, less than max_cap.
#' @param max_cap = The maximum capture rate to consider, ranging 0-1, greater than min_cap.
#' @param min_stor = The minimum storage rate to consider, ranging 0-1, less than max_stor.
#' @param max_stor = The maximum storage rate to consider, ranging 0-1, greater than min_stor.
#' @param min_turn = The minimum turnover period to consider, in years, less than max_turn.
#' @param max_turn = The maximum turnover period to consider, in years, greater than min_turn.
#' @param fines = Boolean indicating whether to simulate fines or gravels deposits.
#' @param storage_index is the age value index for storage
#' @param source_index is the age value index corresponding to the PMF source_prob
#' @param source_prob is the PMF of debris-flow deposit ages
#' @param gravel_index is the age value index corresponding to the PMF gravel_prob
#' @param gravel_prob is the PMF of gravel ages
#' @param obs The vector of observed deposit ages.
#' @param fines Boolean indicating whether the observed deposits are fines.
#' @return table of mean hit % of test thresholds, using test_fit()
#' @export
#' @seealso fines_synth_n
#' @seealso gravel_synth_n
#' @seealso test_fit
fluvial_fit <- function(batch = 10, n = 10,
                       min_cap = 0, max_cap = 1,
                       min_stor = 0, max_stor = 1,
                       min_turn = 50, max_turn = 10000,
                       storage_index, source_index, source_prob,
                       gravel_index, gravel_prob,
                       obs, fines = F) {
  capture_rates <- runif(batch, min_cap, max_cap)
  storage_rates <- runif(batch, min_stor, max_stor)
  turnovers <- runif(batch, min_turn, max_turn)
  mat <- 0
  if (fines) {
    mat <- t(parallel::mcmapply(function(a,b,c,d,e,f,g,h,i,j) fines_fit_n(a, b, c, d, e, f, g, h, i, j),
                                b = capture_rates, c = storage_rates, d = turnovers,
                                MoreArgs = list(a = n, e = storage_index, f = source_index,
                                                g = source_prob, h = gravel_index,
                                                i = gravel_prob, j = obs)))
  } else {
    mat <- t(parallel::mcmapply(function(a,b,c,d,e,f,g,h) gravel_fit_n(a, b, c, d, e, f, g, h),
                                b = capture_rates, c = storage_rates, d = turnovers,
                                MoreArgs = list(a = n, e = storage_index, f = source_index,
                                                g = source_prob, h = obs)))
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

