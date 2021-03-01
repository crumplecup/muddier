
#' accumulate
#'
#' bootstrap wrapper for recorder function
#'
#' @param n is the number of accumulated arrivals
#' @param ri is the rate in, a float
#' @param ro is the rate out, a float
#' @param ia is a vector of inherited ages, integers
#' @param it is the number of iterations of the bootstrap, an integer
#' @return the median deposit ages of `n` accumulations after `it` bootstraps
#' @export

accumulate <- function(n, ri, ro, ia, it = 100) {
  ar <- array(0, c(n, it))
  for (i in 1:it) {
    b <- recorder(n, ri, ro, ia)
    ar[1:length(b), i] <- b
  }
  vec <- apply(ar, 1, median)
  return(vec)
}


#' accumulate1
#'
#' bootstrap wrapper for recorder function
#' adds minimum age
#'
#' @param n is the number of accumulated arrivals
#' @param ri is the rate in, a float
#' @param ro is the rate out, a float
#' @param ia is a vector of inherited ages, integers
#' @param mn is the minimum age
#' @param it is the number of iterations of the bootstrap, an integer
#' @return the median deposit ages of `n` accumulations after `it` bootstraps
#' @seealso recorder1
#' @export

accumulate1 <- function(n, ri, ro, ia, mn, it = 100) {
  ar <- array(0, c(n, it))
  for (i in 1:it) {
    b <- recorder1(n, ri, ro, ia, mn)
    ar[1:length(b), i] <- b
  }
  vec <- apply(ar, 1, median)
  return(vec)
}


#' accumulater
#'
#' bootstrap wrapper for accumulate function
#'
#' @param n is the number of accumulated arrivals
#' @param ti_rng is the range of rates in, a vector of floats
#' @param to_rng is the range of rates out, a vector of floats
#' @param ia is a vector of inherited ages, integers
#' @param obs is the vector of observed estimated mean charcoal ages
#' @param batch is the batch size of random draws
#' @param it is the number of iterations of the bootstrap, an integer
#' @return the median deposit ages of `n` accumulations after `it` bootstraps
#' @seealso accumulate
#' @seealso recorder
#' @export

accumulater <- function(n, ti_rng, to_rng, ia, obs, batch = 10, it = 100) {
  go <- runif(batch, min(to_rng), max(to_rng))
  gi <- runif(batch, min(ti_rng), max(ti_rng))

  res <- mcmapply(
    function(x, y, a, b, c) accumulate(a, x, y, b, c),
    x = gi, y = go,
    MoreArgs = list(a = n, b = ia, c = it)
  )
  gof <- apply(res, 2, function(x) cdf_gof(x, obs))
  gof <- matrix(unlist(gof), ncol = 2, byrow = TRUE)
  df <- data.frame(
    ti = gi,
    to = go,
    ks = gof[, 1],
    kp = gof[, 2]
  )
  return(df)
}


#' accumulater1
#'
#' bootstrap wrapper for accumulate function
#' adds minimum age
#'
#' @param n is the number of accumulated arrivals
#' @param ti_rng is the range of rates in, a vector of floats
#' @param to_rng is the range of rates out, a vector of floats
#' @param ia is a vector of inherited ages, integers
#' @param obs is the vector of observed estimated mean charcoal ages
#' @param batch is the batch size of random draws
#' @param it is the number of iterations of the bootstrap, an integer
#' @return the median deposit ages of `n` accumulations after `it` bootstraps
#' @seealso accumulate
#' @seealso recorder
#' @export

accumulater1 <- function(n, ti_rng, to_rng, ia, obs, batch = 10, it = 100) {
  go <- runif(batch, min(to_rng), max(to_rng))
  gi <- runif(batch, min(ti_rng), max(ti_rng))
  mn <- min(obs)

  res <- mcmapply(
    function(x, y, a, b, c, d) accumulate1(a, x, y, b, c, d),
    x = gi, y = go,
    MoreArgs = list(a = n, b = ia, c = mn, d = it)
  )
  gof <- apply(res, 2, function(x) cdf_gof(x, obs))
  gof <- matrix(unlist(gof), ncol = 2, byrow = TRUE)
  df <- data.frame(
    ti = gi,
    to = go,
    ks = gof[, 1],
    kp = gof[, 2]
  )
  return(df)
}


#' boot_accum
#'
#' bootstrap utility for estimating confidence intervals
#' wraps accumulater()
#'
#' @param n is the number of accumulated arrivals
#' @param ti_rng is the range of rates in, a vector of floats
#' @param to_rng is the range of rates out, a vector of floats
#' @param ia is a vector of inherited ages, integers
#' @param obs is the vector of observed estimated mean charcoal ages
#' @param batch is the batch size of random draws
#' @param it is the number of iterations of the bootstrap, an integer
#' @return the result with the lowest ks test value
#' @seealso accumulater
#' @seealso accumulate
#' @seealso recorder
#' @export

boot_accum <- function(n, ti_rng, to_rng, ia, obs, batch = 1000, it = 100) {
  robs <- sample(obs, length(obs), replace = TRUE)
  rec <- accumulater(n, ti_rng, to_rng, ia, robs, batch, it)
  res <- rec[rec$ks == min(rec$ks), ]
  res[1, ]
}



#' recorder
#'
#' models arrivals and removals from a reservoir
#'
#' @param n is the number of accumulated arrivals
#' @param ri is the rate in, a float
#' @param ro is the rate out, a float
#' @param ia is a vector of inherited ages, integers
#' @return deposit ages of accumulated arrivals, a numeric vector
#' @export

recorder <- function(n, ri, ro, ia)  {
  im <- 0
  k <- 1
  rec <- 0
  om <- 0
  too_old <- 50000

  while (length(rec) < n &
         im < too_old
  )  {
    om <- om + rexp(1, ro)
    while (
      im < om &
      length(rec) <= n
    ) {
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
  #  rec <- rec[1:k]
  rec <- sort(max(rec) - rec)
  rc <- 0
  for (i in seq_along(rec)) {
    rc[i] <- rec[i] + sample(ia, 1)
  }
  return(sort(rc - min(rc)))
}



#' recorder1
#'
#' models arrivals and removals from a reservoir
#' adds a minimum age
#'
#' @param n is the number of accumulated arrivals
#' @param ri is the rate in, a float
#' @param ro is the rate out, a float
#' @param ia is a vector of inherited ages, integers
#' @param mn is a minimum age
#' @return deposit ages of accumulated arrivals, a numeric vector
#' @export

recorder1 <- function(n, ri, ro, ia, mn)  {
  im <- 0
  k <- 1
  rec <- 0
  om <- 0
  too_old <- 50000

  while (length(rec) < n &
         im < too_old
  )  {
    om <- om + rexp(1, ro)
    while (
      im < om &
      length(rec) <= n
    ) {
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
  #  rec <- rec[1:k]
  rec <- sort(max(rec) - rec)
  rc <- 0
  for (i in seq_along(rec)) {
    rc[i] <- rec[i] + sample(ia, 1)
  }
  return(sort(rc - min(rc) + mn))
}


#' rec_vol
#'
#' function to record volumes
#'
#' @param node is a spatial object (creek node)
#' @param n is an integer representing the number of years to accumulate
#' @param vi is the volume input rate, a float
#' @param vo is the volume output rate, a float
#' @param return matrix of volumes and levels
#' @export

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


#' fit_volumes
#'
#' function to fit volumes to creek nodes
#'
#' @param nodes is a spatial object (creek nodes)
#' @param n is an integer representing years to simulate accumulation
#' @param vi is the volumetric input rate, a float
#' @param vo is the volumetric output rate, a float
#' @param it is an integer representing the number of iterations per node
#' @return a data.frame with cols c('vol', 'lvl')
#' @export

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




#' select_rate_pairs
#'
#' function to return rate pairs given a preference matrix
#'
#' @param ar is the preference matrix (a pdf superimposed on a grid)
#' @param xs is the vector of values along the x-axis of `ar`
#' @param ys is the vector of values along the y-axis of `ar`
#' @return a pair of floats c(x, y)
#' @export
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


#' weight_array
#'
#' Returns preference matrix given topology matrix.
#' The scale parameter deflates preference weighting from zero to one,
#' inflates preference weighting above one,
#' and crashes below zero.
#'
#' @param ar is a topology matrix from weight_grid()
#' @param scl is the scale parameter
#' @return a preference matrix
#' @seealso weight_grid
#' @export

weight_array <- function(ar, scl = 1) {
  vec <- as.vector(ar)
  wmax <- max(vec)
  p <- (1 - vec) * scl * (1 / length(vec))
  p <- p / sum(p)
  array(p, dim(ar))
}


#' weight_grid
#'
#' Returns topology of coordinate on grid given record of test values.
#' Types for `test`:
#' Kolmogorov-Smirnov statistic = 'ks'
#' Kuiper statistic = 'kp'
#' `xy_grid` determines the size of the outgoing topology array.
#' Types for `by`:
#' Unit cross-sectional area in m2 = 'vol'
#' Average cross-sectional depth per node = 'lvl'
#' `scl` ranges [0,1] and deflates the effects of neighbors on each other.
#'
#' @param pos is the coordinate pair of data in the record
#' @param xy_grid is the topology array to replace
#' @param rec is the record of test values
#' @param test is a string naming the test type c('ks', 'kp)
#' @param by is a string naming the test varaible c('vol', 'lvl')
#' @param scl is a scaling parameter, a float
#' @return the estimated test value at `pos`
#' @export

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


#' update_search_array
#'
#' Returns preference topology given a record of test values.
#' Types for `test`:
#' Kolmogorov-Smirnov statistic = 'ks'
#' Kuiper statistic = 'kp'
#' Types for `by`:
#' Unit cross-sectional area in m2 = 'vol'
#' Average cross-sectional depth per node = 'lvl'
#'
#' @param rec is a data.frame of test values with cols c(x, y, ks, kp)
#' @param xs is a vector of values for the x-axis of topology array
#' @param ys is a vector of values for the y-axis of topology array
#' @param ar is a preference topology to update
#' @param test is a string naming the test type c('ks', 'kp)
#' @param by is a string naming the test varaible c('vol', 'lvl')
#' @param d_scl is the scale value for weight_grid()
#' @param w_scl is the scale value for weight_array()
#' @return a preference topology array

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


#' search_topology
#'
#'
#'
#' @param so is a spatial object (creek nodes)
#' @param x_range is a vector of values for the x-axis of the preference array
#' @param y_range is a vector of values for the y-axis of the preference array
#' @param rec is a data.frame of test values with cols c(x, y, ks, kp)
#' @param bt is the number of times to simulate accumulation over `so`
#' @param it is the number of times to simulate accumulation per node of `so`
#' @param n is the number of years to accumulate deposits
#' @param test is a string naming the test type c('ks', 'kp)
#' @param by is a string naming the test varaible c('vol', 'lvl')
#' @param grid is the number of intervals to divide x and y ranges into
#' @param d_scl is the scale value for weight_grid()
#' @param w_scl is the scale value for weight_array()
#' @return a record of test value results at randomly generated coords in `x_range` and `y_range`
#' @export

search_topology <- function(so,
                            x_range,
                            y_range,
                            rec = 0,
                            bt = 10,
                            it = 20,
                            n = 10000,
                            test = 'ks',
                            by = 'vol',
                            grid = 100,
                            d_scl = 1,
                            w_scl = 3,
                            linked = FALSE,
                            backfill = TRUE
                            ) {
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
    out <- 0
    if (length(vec) > 1) {
      ar <- update_search_array(
        matrix(unlist(vec), ncol = 7, byrow = TRUE),
        xs, ys, ar, test, by, d_scl, w_scl)
      present_topology(ar, xs, ys)
    }
    rates <- select_rate_pairs(ar, xs, ys)
    if (linked) {
      res <- fit_volumes1(so, n, rates[1], rates[2], it)
      boot <- res[[1]]
      out <- res[[2]]
    }
    if (backfill) {
      res <- backfill(so, n, rates[1], rates[2], it)
      boot <- res[[1]]
      out <- res[[2]]
    }
    if (!linked & !backfill) boot <- fit_volumes(so, n, rates[1], rates[2], it)
    vol_gof <- cdf_gof(boot[, 1], so$xsec_area)
    lvl_gof <- cdf_gof(boot[, 2], so$lvl)

    if (length(vec) > 1) vec <- c(vec, rates, vol_gof, lvl_gof, out)
    if (length(vec) == 1) vec <- c(rates, vol_gof, lvl_gof, out)
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
  matrix(unlist(vec), ncol = 7, byrow = TRUE)
}


#' present_topology
#'
#' prints a visualization of the topology array
#'
#' @param ar is a topology array
#' @param xs is a vector of values for the x-axis of the array
#' @param ys is a vector of values for the y-axis of the array
#' @return prints a visualization of the topology array
#' @export

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


#' weight_by_delprob
#'
#' converts MB delivery probabilities to Rose-Lancaster optimized delivery probabilites
#'
#' @param delprob is a vector of delivery probabilities, f64
#' @param weight is a matrix with cols c(cdf, delprob)
#' @return vector of optimized delivery probabilities
#' @export

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


#' record_volume
#'
#' Record acculumation at a node.  Takes arrivals and emits departures.
#' Nodes are linked, with deposits departing downstream.
#'
#' @param node is a row from a spatial dataframe (creek nodes)
#' @param n is the number of years to simulate accumulation
#' @param vi is the volumetric input rate
#' @param vo is the volumetric output rate
#' @param arrivals is a dataframe describing upstream arrivals
#' @return a list with elements c(accumulation, departures)
#' @seealso fit_volumes1
#' @export

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

#' fit_volumes1
#'
#' Function to fit volume rate to observed volumes.
#' Creek nodes are linked, passing deposits downstream.
#'
#' @param nodes is a spatial object (creek nodes)
#' @param n is the number of years to simulate accumulation, an integer
#' @param vi is the volumetric input rate
#' @param vo is the volumetric output rate
#' @param it is the number of iterations to simulate accumulation, an integer
#' @return the median accumulation record after `it` iterations in volume and level
#' @seealso record_volume
#' @export

fit_volumes1 <- function(nodes, n, vi, vo, it = 10) {
  vol <- array(0, c(nrow(nodes), it))
  lvl <- array(0, c(nrow(nodes), it))
  out <- 0
  nodes <- nodes[order(nodes$ToMouth_km, decreasing = TRUE), ]
  arrivals <- data.frame(t = 0, r = 'input', vol = 0, lvl = 0)
  res <- mclapply(1:it, function(x) voluminous(nodes, n, vi, vo, arrivals))
  for (i in 1:it) {
    rec <- res[[i]]
    vol[ , i] <- rec[[1]]
    lvl[ , i] <- rec[[2]]
    dep <- rec[[3]]
    dep <- dep[dep$t < n, ]
    out[i] <- sum(dep$vol) / n
  }

  df <- data.frame(
    vol = apply(vol, 1, median),
    lvl = apply(lvl, 1, median)
    )
  return(list(df, mean(out)))
}



#' voluminous
#'
#' Wrapper for record_volume.
#'
#' @param nodes is a spatial object (creek nodes)
#' @param n is the number of years to simulate accumulation, an integer
#' @param vi is the volumetric input rate
#' @param vo is the volumetric output rate
#' @param arrivals is the record of arrivals, a data.frame with cols (t, r, vol, lvl)
#' @return a list with elements c(volumes, levels, arrivals)
#' @export

voluminous <- function(nodes, n, vi, vo, arrivals) {
  vol <- 0
  lvl <- 0
  for (i in 1:nrow(nodes)) {
    res <- record_volume(nodes[i, ], n, vi, vo, arrivals)
    rec <- res[[1]]
    arrivals <- res[[2]]
    vol[i] <- rec[1,1]
    lvl[i] <- rec[1,2]
  }
  return(list(vol = vol, lvl = lvl, arrivals = arrivals))
}



#' voluminous1
#'
#' Simulates accumulation record, using a linked-bucket model
#' with backfilling based on average elevation of bank deposits.
#'
#' @param nodes is a spatial object (creek nodes)
#' @param n is the number of years to simulate accumulation, an integer
#' @param vi is the volumetric input rate
#' @param vo is the volumetric output rate
#' @param b_scl is the proportion of volume that settles out of deposits at obstructions
#' @return a list with elements c(volumes, levels, arrivals)
#' @seealso eventer
#' @export

voluminous1 <- function(nodes, n, vi, vo, b_scl) {
  nodes <- nodes[order(nodes$ToMouth_km), ]
  vol <- array(0, nrow(nodes))
  lvl <- array(0, nrow(nodes))
  elev <- nodes$elev
  arr <- mapply(function(x, y, z, n, m) as.matrix(eventer(x, y, z, n, m)),
                x = nodes$ri, y = nodes$NODE_ID,
                MoreArgs = list(z = vi, n = n, m = 0))
  rec <- arr[[1]]
  for (i in 2:length(arr)) {
    rec <- rbind(rec, arr[[i]])
  }
  dep <- mapply(function(x, y, z, n, m) as.matrix(eventer(x, y, z, n, m)),
                x = nodes$ro, y = nodes$NODE_ID,
                MoreArgs = list(z = vo, n = n, m = 1))
  for (i in seq_along(dep)) {
    rec <- rbind(rec, dep[[i]])
  }
  rec <- rec[order(rec[, 1]), ]
#  print('record:')
#  print(head(rec[rec[,1] > 0, ], 50))
  ids <- nodes$NODE_ID
  departures <- data.frame(age = 0, vol = 0, type = 0, id = 0)
  for (i in 1:nrow(rec)) {
#    print(paste0('on record ', i))
    k <- sum((ids == rec[i, 4]) * 1:length(ids))
#    print(paste0('on index ', k))
#    print(paste0('is input: ', rec[i, 3] == 0))
    if (rec[i, 3] == 0) {
      vol[k] <- vol[k] + rec[i, 2]
      lev <- lvl[k] + rec[i, 2] /
        as.numeric(nodes$valley_width[k])
      lvl[k] <- lev
      elev[k] <- elev[k] + lev
    }
#    print(paste0('is output: ', rec[i, 3] == 1))
    if (rec[i, 3] == 1) {
      scoop <- 0
      settle <- 0
      if (rec[i, 2] > vol[k]) {
        scoop <- vol[k]
      }
      if (rec[i, 2] <= vol[k]) {
        scoop <- rec[i, 2]
      }
#      print(paste0('scoop is ', scoop))
      og_k <- k
      vol[k] <- vol[k] - scoop
      lvl[k] <- lvl[k] - scoop / as.numeric(nodes$valley_width[k])
      elev[k] <- elev[k] - scoop / as.numeric(nodes$valley_width[k])

      going <- TRUE
      idl <- length(ids)

      while (k <= idl & going == TRUE) {
        if (k == idl) {
          rec[i, 2] <- scoop
          departures <- rbind(departures, rec[i, ])
          going <- FALSE
#          print('off and away')
        }
        if (k < idl) {
          k <- k + 1
#          print(paste0('going to ', k))
          if (lvl[k-1] < lvl[k]) {
            settle <- scoop * b_scl
            scoop <- scoop - settle
            vol[k] <- vol[k] + settle
            lvl[k] <- lvl[k] + settle / as.numeric(nodes$valley_width[k])
            elev[k] <- elev[k] + settle / as.numeric(nodes$valley_width[k])
#            print('never got started')
          }
          if (lvl[k-1] >= lvl[k]) {
            r <- runif(1)
            if (r <= nodes$rdp[k]) {
              settle <- scoop * b_scl
              scoop <- scoop - settle
              vol[k] <- vol[k] + settle
              lvl[k] <- lvl[k] + settle / as.numeric(nodes$valley_width[k])
              elev[k] <- elev[k] + settle / as.numeric(nodes$valley_width[k])
#              print(paste0('new vol is: ', vol[k]))
            }
          }
        }
      }
    }
  }
  return(list(vol = vol, lvl = lvl, elev = elev, departures = departures)
  )
}



#' eventer
#'
#' Utility for voluminous().
#' Generates an event data.frame with cols (age, vol, type, id),
#' given a temporal and volumetric event rate, over a given period of time.
#'
#' @param rt is the temporal event rate
#' @param id is the creek node number, an integer.
#' @param vl is the volumetric event rate
#' @param n is the period of time to simulate events
#' @param lab is the type label, which sounds like a string but is numeric!
#' @return a data.frame of events with cols (age, vol, type, id)
#' @export

eventer <- function(rt, id, vl, n, lab) {
  t <- 0
  age <- 0
  vol <- 0
  k <- 1
  while(t < n) {
    t <- t + rexp(1, rt)
    if (t <= n) {
      age[k] <- t
      vol[k] <- rexp(1, vl)
      k <- k + 1
    }
  }
  return(
    data.frame(
      age = age,
      vol = vol,
      type = rep(lab, length(vol)),
      id = rep(id, length(vol))
    )
  )
}


#' backfill
#'
#' Function to fit volume rate to observed volumes.
#' Creek nodes are linked, passing deposits downstream.
#' Deposits backfill based upon average elevation
#'
#' @param nodes is a spatial object (creek nodes)
#' @param n is the number of years to simulate accumulation, an integer
#' @param vi is the volumetric input rate
#' @param vo is the volumetric output rate
#' @param it is the number of iterations to simulate accumulation, an integer
#' @param b_scl is the proportion of volume that settles out of deposits at obstructions
#' @return the median accumulation record after `it` iterations in volume and level
#' @seealso record_volume
#' @export

backfill <- function(nodes, n, vi, vo, b_scl, it = 10) {
  vol <- array(0, c(nrow(nodes), it))
  lvl <- array(0, c(nrow(nodes), it))
  elev <- array(0, c(nrow(nodes), it))
  out <- 0
  nodes <- nodes[order(nodes$ToMouth_km, decreasing = TRUE), ]
  arrivals <- data.frame(t = 0, r = 'input', vol = 0, lvl = 0)
  res <- lapply(1:it, function(x) voluminous1(nodes, n, vi, vo, b_scl))
  for (i in 1:it) {
    rec <- res[[i]]
    vol[ , i] <- rec[[1]]
    lvl[ , i] <- rec[[2]]
    elev[ , i] <- rec[[3]]
    dep <- rec[[4]]
    dep <- dep[dep[,1] < n, ]
    out[i] <- sum(dep$vol) / n
  }

  df <- data.frame(
    vol = apply(vol, 1, median),
    lvl = apply(lvl, 1, median),
    elev = apply(elev, 1, median)
  )
  return(list(record = df, output = mean(out)))
}


#' slope_to_elev
#'
#' Converts slopes to elevations based on distance from a base.
#'
#' @param slope is a numeric vector of slopes at creek nodes
#' @param dist is a numeric vector of distances to mouth at creek nodes
#' @param base is the starting elevation
#' @return a vector of elevations
#' @export

slope_to_elev <- function(slope, dist, base = 0) {
  elev <- base
  for (i in seq_along(slope)) {
    if (i == 1) {
      elev[i] <- base
    }
    if (i > 1) {
      elev[i] <- elev[i-1] + ((dist[i] - dist[i-1]) * slope[i-1] * 1000)
    }
  }
  elev
}


#' n_dimensional_search
#'
#' For finding a needle in a haystack.
#' Uses backfill() to simulate accumulation.
#'
#' @param so is a spatial object (creek nodes)
#' @param ti_rng is the range of temporal input rates to search
#' @param to_rng is the range of temporal output rates to search
#' @param vi_rng is the range of volumetric input rates to search
#' @param vo_rng is the range of volumetric input rates to search
#' @param ti_scl is the range of temporal input scales to search
#' @param to_scl is the range of temporal output scales to search
#' @param ks_scl is the range of streampower scales to search
#' @param dp_scl is the range of deposit probability scales to search
#' @param b_scl is the proportion of volume that settles out of deposits at obstructions
#' @param it is the number of iterations to perform at each point
#' @param n is the number of years to simulate accumulation, an integer
#' @seealso backfill
#' @return a list with elements (record, mean output rate)
#' @export

n_dimensional_search <- function(so,
                                 ti_rng,
                                 to_rng,
                                 vi_rng,
                                 vo_rng,
                                 ti_scl = 1,
                                 to_scl = 1,
                                 ks_scl = 1,
                                 dp_scl = 1,
                                 m_scl = 1,
                                 b_scl = 1,
                                 type = 0,
                                 it = 100,
                                 n = 10000,
                                 batch = 10) {
  go <- runif(batch, min(to_rng), max(to_rng))
  gi <- runif(batch, go, max(ti_rng))
  vo <- runif(batch, min(vo_rng), max(vo_rng))
  vi <- runif(batch, min(vi_rng), max(vi_rng))
  to_scale <- runif(batch, min(to_scl), max(to_scl))
  ks_scale <- runif(batch, min(ks_scl), max(ks_scl))
  ti_scale <- runif(batch, min(ti_scl), max(ti_scl))
  dp_scale <- runif(batch, min(dp_scl), max(dp_scl))
  m_scale <- runif(batch, min(m_scl), max(m_scl))
  b_scale <- runif(batch, min(b_scl), max(b_scl))
  if (type == 0) type <- sample(c('ks', 'dp', 'mixed'), 1)
  res <- mcmapply(
    function (x, a, b, c, d, e, f, g, h, j, k, l, m, p) backfill(
      rater(x, a, b, c, d, e, f, g, h), l, j, k, m, p),
    a = gi, b = go, c = ti_scale, d = to_scale, e = ks_scale, f = dp_scale, g = m_scale, j = vi, k = vo, m = b_scale,
    MoreArgs = list(x = so, h = type, l = n, p = it)
  )


  vol_gof <- array(0, c(batch, 2))
  lvl_gof <- array(0, c(batch, 2))
  elev_gof <- array(0, c(batch, 2))
  out <- array(0, batch)
  vor <- array(0, batch)
  for (i in 1:(length(res)/2)) {
    bt <- res[[i + (i-1)]]
    vol_gof[i, ] <- as.numeric(cdf_gof(as.numeric(unlist(bt[1])), so$xsec_area))
    lvl_gof[i, ] <- as.numeric(cdf_gof(as.numeric(unlist(bt[2])), so$xsec_area / as.numeric(so$valley_width)))
    elev_gof[i, ] <- as.numeric(cdf_gof(as.numeric(unlist(bt[3])), slope_to_elev(so$slope, so$ToMouth_km)))
    out[i] <- res[[i + (i-1) + 1]]
  }
  data.frame(
    ti = gi,
    to = go,
    vi = vi,
    vo = vo,
    ti_scl = ti_scale,
    to_scl = to_scale,
    ks_scl = ks_scale,
    dp_scl = dp_scale,
    m_scl = m_scale,
    b_scl = b_scale,
    ks_vol = vol_gof[ ,1],
    kp_vol = vol_gof[ ,2],
    ks_lvl = lvl_gof[ ,1],
    kp_lvl = lvl_gof[ ,2],
    ks_elev = elev_gof[, 1],
    kp_elev = elev_gof[ ,2],
    out = out,
    type = type,
    stringsAsFactors = FALSE
  )
}


#' rater
#'
#' Derives variables related to temporal input and output rates.
#' The average rate per node is the global rate divided by the number of sampled nodes in sa_nodes.
#' Types of `type` include c('ks', 'dp', 'mixed'):
#' Streampower coefficient = 'ks'.
#' Inverse delivery probability = 'dp'.
#' Both = 'mixed'.
#'
#' @param so is a spatial object (creek nodes)
#' @param gi is a global temporal input rate
#' @param go is a global temporal output rate
#' @param ti_scl is the range of temporal input scales to search
#' @param to_scl is the range of temporal output delivery scales to search
#' @param ks_scl is the range of temporal output streampower scales to search
#' @param dp_scl is the range of deposit probability scales to search
#' @param m_scl is the range of distance probability scales to search
#' @param type is a string indicated the type of output weight to apply
#' @param sa_nodes is the spatial object of sampled nodes
#' @return spatial object with rates calculated
#' @seealso n_dimensional_search
#' @export

rater <- function(so, gi, go, ti_scl, to_scl, ks_scl, dp_scl, m_scl, type, sa = sa_nodes) {
  so <- so[order(so$ToMouth_km, decreasing = FALSE), ]
  wt <- unlist(lapply(so$DebrisFlow, function(x) weight_by_delprob(x)))
  dp <- so$DebrisFlow * wt
  rdp <- dp * dp_scl
  ldp_k <- exp(log(dp) - mean(log(dp)))
#  wdp <- log(dp) - min(log(dp))
  lvl <- so$xsec_area / as.numeric(so$valley_width)
  elev <- slope_to_elev(so$slope, so$ToMouth_km)
  pwr_df <- data.frame(CA = log(so$contr_area),
                       slope = log(so$slope))
  pwr_mod <- lm(CA ~ slope, data = pwr_df)
  pwr_prd <- predict(pwr_mod, newdata = data.frame(slope = log(so$slope)))
  pwr_k <- log(so$contr_area) - pwr_prd
  wdp <- -pwr_k
  wdp <- wdp - min(wdp)
  mth <- so$ToMouth_km - min(so$ToMouth_km)
  mth <- mth / max(mth)
  mth <- 1 - (mth - m_scl * mth^2)

  to_ave <- go / nrow(sa)
  ti_ave <- gi / nrow(sa)
  to <- 0
  if (type == 'ks') to <- to_ave * exp(pwr_k) * ks_scl
  if (type == 'dp') to <- to_ave * (1 / ldp_k) * to_scl
  if (type == 'mixed') to <- to_ave * exp(pwr_k) * (1 / ldp_k) * ks_scl * to_scl
  to[to <= 0] <- min(to[to > 0])
  to <- to * sum(to_ave * nrow(so)) / sum(to)
  ti <- ti_ave * wdp * mth * ti_scl
  ti[ti <= 0] <- min(ti[ti > 0])
  ti <- ti * sum(ti_ave * nrow(so)) / sum(ti)
  so$dp <- dp
  so$rdp <- rdp
  so$ri <- ti
  so$ro <- to
  so$lvl <- lvl
  so$elev <- elev
  so
}



#' voluminous2
#'
#' Simulates accumulation record, using a linked-bucket model
#' with backfilling based on average elevation of bank deposits.
#'
#' @param nodes is a spatial object (creek nodes)
#' @param n is the number of years to simulate accumulation, an integer
#' @param vi is the volumetric input rate
#' @param vo is the volumetric output rate
#' @param tp is the turbulent deposition probability
#' @param lp is the laminar deposition probability
#' @return a list with elements c(volumes, levels, arrivals)
#' @seealso eventer
#' @export

voluminous2 <- function(nodes, n, vi, vo, tp = 0.1, lp = 0.5) {
  nodes <- nodes[order(nodes$ToMouth_km, decreasing = TRUE), ]
  vol <- array(0, nrow(nodes)) # unit cross-sectional area m2
  lvl <- array(0, nrow(nodes)) # average valley depth m
  lam <- array(0, nrow(nodes)) # laminar or turbulent boolean
  elev <- nodes$elev

  recs <- 1:nrow(nodes)
  recs <- lapply(recs, function(x) data.frame(id = x, vol = 0, age = 0, type = 0))

  # assemble record of arrivals and departures
  arr <- mapply(function(x, y, z, n, m) as.matrix(eventer(x, y, z, n, m)),
                x = nodes$ri, y = nodes$NODE_ID,
                MoreArgs = list(z = vi, n = n, m = 0))
  rec <- arr[[1]]
  for (i in 2:length(arr)) {
    rec <- rbind(rec, arr[[i]])
  }
  dep <- mapply(function(x, y, z, n, m) as.matrix(eventer(x, y, z, n, m)),
                x = nodes$ro, y = nodes$NODE_ID,
                MoreArgs = list(z = vo, n = n, m = 1))
  for (i in seq_along(dep)) {
    rec <- rbind(rec, dep[[i]])
  }
  rec <- rec[order(rec[, 1]), ]
  ids <- nodes$NODE_ID
  departures <- data.frame(id =0, vol = 0, age = 0, type = 0)

  # process arrivals and departures
  for (i in 1:nrow(rec)) {
    # k is id in ids of node on which rec occurs
    k <- sum((ids == rec[i, 4]) * 1:length(ids))
    # deposit input
    if (rec[i, 3] == 0) {
      vol[k] <- vol[k] + rec[i, 2]
      lev <- lvl[k] + rec[i, 2] /
        as.numeric(nodes$valley_width[k])
      lvl[k] <- lev
      recs[[k]] <- rbind(recs[[k]], c(k, rec[i, 2], rec[i, 1], 0))
      # mark upstream portions laminar (inundated)
      # if their elevation is below the height of the deposit
      thresh <- elev[k] + lvl[k]
      l <- k
      if (l > 1) l <- l - 1
      while ((elev[l] + lvl[l]) <= thresh & l > 1) {
        lam[l] <- 1
        l <- l - 1
      }
    }
    # withdraw output
    if (rec[i, 3] == 1 & vol[k] > 0) {
      scoop <- 0
      if (rec[i, 2] > vol[k]) {
        scoop <- vol[k]
      }
      if (rec[i, 2] <= vol[k]) {
        scoop <- rec[i, 2]
      }
      vol[k] <- vol[k] - scoop
      lvl[k] <- lvl[k] - scoop / as.numeric(nodes$valley_width[k])
      going <- TRUE
      idl <- length(ids)

      # update node volume record
      slug <- data.frame(id = k, vol = 0, age = 0, type = 1)
      cup <- 0
      nrec <- recs[[k]]
      while (sum(slug[, 2]) < scoop) {
        strat <- sample(1:nrow(nrec), 1)
        cup <- scoop - sum(slug[, 2])
        if (cup <= nrec[strat, 2]) { # scoop takes portion of strat
          slug <- rbind(slug, c(k, cup, nrec[strat, 3], 1))
          nrec[strat, 2] <- nrec[strat, 2] - cup
        }
        if (cup > nrec[strat, 2]) { # scoop consumes strat
          slug <- rbind(slug, c(k, nrec[strat, 2], nrec[strat, 3], 1))
        }
      }
      recs[[k]] <- nrec
      slug <- slug[slug[, 2] > 0, ] # drop zero volumes from slug

      # send output downstream
      while (k <= idl & going == TRUE) {
        # case: end of study area
        if (k == idl) {
          departures <- rbind(departures, slug)
          going <- FALSE
        }
        # case: going downstream
        if (k < idl) {
          k <- k + 1
          r <- runif(1)
          settle <- r <= tp # turbulent probability
          if (lam[k]) settle <- r <= lp # laminar probability
          if (settle) {
            lump <- sum(slug[, 2])
            vol[k] <- vol[k] + lump
            lvl[k] <- lvl[k] + lump / as.numeric(nodes$valley_width[k])
            recs[[k]] <- rbind(recs[[k]], slug)
            going <- FALSE
          }
        }
      }
    }
  }
  rex <- 0
  for (i in seq_along(recs)) {
    if (i == 1) rex <- recs[[i]]
    if (i > 1 ) rex <- rbind(rex, recs[[i]])
  }
  return(list(vol = vol, lvl = lvl, recs = rex, dep = departures)
  )
}

#' backfill1
#'
#' Function to fit volume rate to observed volumes.
#' Creek nodes are linked, passing deposits downstream.
#' Deposits backfill based upon average elevation
#'
#' @param nodes is a spatial object (creek nodes)
#' @param n is the number of years to simulate accumulation, an integer
#' @param vi is the volumetric input rate
#' @param vo is the volumetric output rate
#' @param tp is the turbulent deposition probability
#' @param lp is the laminar deposition probability
#' @param it is the number of iterations to simulate accumulation, an integer
#' @return the median accumulation record after `it` iterations in volume and level
#' @seealso record_volume
#' @export

backfill1 <- function(nodes, n, vi, vo, tp, lp, it = 10) {
  vol <- array(0, c(nrow(nodes), it))
  lvl <- array(0, c(nrow(nodes), it))
  out <- 0
  nodes <- nodes[order(nodes$ToMouth_km, decreasing = TRUE), ]
  arrivals <- data.frame(t = 0, r = 'input', vol = 0, lvl = 0)
  res <- lapply(1:it, function(x) voluminous2(nodes, n, vi, vo, tp, lp))
  for (i in 1:it) {
    rec <- res[[i]]
    vol[ , i] <- rec[[1]]
    lvl[ , i] <- rec[[2]]
    recs <- rec[[3]]
    df_vol <- sum(recs$vol[recs$type == 0])
    df_pct <- df_vol / sum(recs$vol)
    dep <- rec[[4]]
    dep <- dep[dep[,3] < n, ]
    out[i] <- sum(dep$vol) / n
  }

  df <- data.frame(
    vol = apply(vol, 1, median),
    lvl = apply(lvl, 1, median)
  )
  return(list(record = df, output = mean(out), df_pct = df_pct))
}



#' ndim_search
#'
#' For finding a needle in a haystack.
#' Uses backfill() to simulate accumulation.
#'
#' @param so is a spatial object (creek nodes)
#' @param ti_rng is the range of temporal input rates to search
#' @param to_rng is the range of temporal output rates to search
#' @param vi_rng is the range of volumetric input rates to search
#' @param vo_rng is the range of volumetric input rates to search
#' @param tp is the turbulent deposition probability
#' @param lp is the laminar deposition probability
#' @param it is the number of iterations to perform at each point
#' @param n is the number of years to simulate accumulation, an integer
#' @seealso backfill
#' @return a list with elements (record, mean output rate)
#' @export

ndim_search <- function(so,
                        ti_rng,
                        to_rng,
                        vi_rng,
                        vo_rng,
                        tp = 0.1,
                        lp = 0.5,
                        it = 100,
                        n = 10000,
                        batch = 10) {
  go <- runif(batch, min(to_rng), max(to_rng))
  gi <- runif(batch, go, max(ti_rng))
  vo <- runif(batch, min(vo_rng), max(vo_rng))
  vi <- runif(batch, min(vi_rng), max(vi_rng))
  tps <- runif(batch, min(tp), max(tp))
  lps <- runif(batch, min(lp), max(lp))
  res <- mcmapply(
    function (x, a, b, c, d, e, f, g, h) backfill1(
      rater1(x, a, b), c, d, e, f, g, h),
    a = gi, b = go, d = vi, e = vo, f = tps, g = lps,
    MoreArgs = list(x = so, c = n, h = it)
  )

  vol_gof <- array(0, c(batch, 2))
  lvl_gof <- array(0, c(batch, 2))
  out <- array(0, batch)
  df_pct <- array(0, batch)
  for (i in 1:(length(res)/3)) {
    bt <- res[[3 * (i-1) + 1]]
    vol_gof[i, ] <- as.numeric(cdf_gof(as.numeric(unlist(bt[1])), so$xsec_area))
    lvl_gof[i, ] <- as.numeric(cdf_gof(as.numeric(unlist(bt[2])), so$xsec_area / as.numeric(so$valley_width)))
    out[i] <- res[[3 * (i-1) + 2]]
    df_pct[i] <- res[[3 * (i-1) + 3]]
  }
  data.frame(
    ti = gi,
    to = go,
    vi = vi,
    vo = vo,
    tp = tps,
    lp = lps,
    ks_vol = vol_gof[ ,1],
    kp_vol = vol_gof[ ,2],
    ks_lvl = lvl_gof[ ,1],
    kp_lvl = lvl_gof[ ,2],
    out = out,
    df_pct = df_pct,
    stringsAsFactors = FALSE
  )
}




#' rater1
#'
#' Derives variables related to temporal input and output rates.
#' The average rate per node is the global rate divided by the number of sampled nodes in sa_nodes.
#'
#' @param so is a spatial object (creek nodes)
#' @param gi is a global temporal input rate
#' @param go is a global temporal output rate
#' @param sa_nodes is the spatial object of sampled nodes
#' @return spatial object with rates calculated
#' @seealso n_dimensional_search
#' @export

rater1 <- function(so, gi, go, sa = sa_nodes) {
  so <- so[order(so$ToMouth_km, decreasing = FALSE), ]
  wt <- unlist(lapply(so$DebrisFlow, function(x) weight_by_delprob(x)))
  dp <- so$DebrisFlow * wt
  dp_wt <- log(dp) - min(log(dp))

  lvl <- so$xsec_area / as.numeric(so$valley_width)
  elev <- slope_to_elev(so$slope, so$ToMouth_km)
  pwr_df <- data.frame(CA = log(so$contr_area),
                       slope = log(so$slope))
  pwr_mod <- lm(CA ~ slope, data = pwr_df)
  pwr_prd <- predict(pwr_mod, newdata = data.frame(slope = log(so$slope)))
  pwr_k <- log(so$contr_area) - pwr_prd
  ca_wt <- so$contr_area / max(so$contr_area)

  to_ave <- go / nrow(sa)
  ti_ave <- gi / nrow(sa)
  to <- to_ave * exp(pwr_k)
  to[to <= 0] <- min(to[to > 0])
  to <- to * sum(to_ave * nrow(so)) / sum(to)

  ti <- ti_ave * dp_wt * ca_wt / exp(pwr_k)
  ti[ti <= 0] <- min(ti[ti > 0])
  ti <- ti * sum(ti_ave * nrow(so)) / sum(ti)
  so$dp <- dp
  so$ri <- ti
  so$ro <- to
  so$lvl <- lvl
  so$elev <- elev
  so
}








