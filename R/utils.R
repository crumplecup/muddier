# Functions for Deriving Inherited Age Distributions


#' boot ids
#'
#' wrapper for lapply, given a rank list from a site, returns
#' `n` bootstraps with replacement of equal length as the rank list
#'
#' @param n is an integer specifying the number of bootstraps
#' @param rl is a rank list output from function `rank_list`
#' @return `n` bootstraps with replacement of length `rl`
#' @seealso rank_list
#' @export

boot_ids <- function(n, rl)  {
  mclapply(rl, function(a) lapply(seq_len(n), function(x,y) sort(sample(y,rep=T)), y = a))
}



#' build_mats condenses a sparse matrix down to a dense one.
#'
#' Given two vectors of equal length, returns a matrix where values of the first vector are nonzero.
#'
#' @param x is a numeric vector of probabilities (a pmf)
#' @param y is a numeric vector of values associated with x
#' @return a dense matrix (x,y) excluding nonzero probabilities in x
#'
#' @export


build_mats <- function(x,y)  {
  mat <- matrix(c(x,y),ncol=2)
  mat <- mat[mat[,1]>0,]
  mat
}


#' derived pmf by subtraction
#'
#' Given two pmfs and an age vector, returns the derived pmf of the convolution y-x.
#' Lengths of x,y and index must be equal.
#'
#' @param x is a numeric vector (younger pmf)
#' @param y is a numeric vector (older pmf)
#' @param index is a numeric vector (years)
#' @return a length(index) numeric vector of the convolved distribution of y-x
#' @seealso convo_minus
#'
#' @export


convo <- function(x,y,index)    {
  probcomb(convo_minus(build_mats(x,index),build_mats(y,index)),sort(index))
}



#' derived pmf by addition
#'
#' Given two pmfs and an age vector, returns the derived pmf of the convolution y+x.
#' Lengths of x,y and index must be equal.
#'
#' @param x is a numeric vector (younger pmf)
#' @param y is a numeric vector (older pmf)
#' @param index is a numeric vector (years)
#' @return a length(index) numeric vector of the convolved distribution of y+x
#' @seealso convo_plus
#'
#' @export


convo_add <- function(x,y,index)    {
  probcomb(convo_plus(build_mats(x,index),build_mats(y,index)),sort(index))
}


#' convolve list
#'
#' returns a list of convolved `pmfs` with ranks in `lis`
#'
#' @param lis is a list of integer ranks
#' @param pmfs is a matrix of pmfs with rownames(vals) (charcoal ages)
#' @return a list of convolved `pmfs` with ranks in `lis`

convo_lis <- function(lis, pmfs = char_pmfs)  {
  index <- as.numeric(rownames(pmfs))
  lapply(lis, function(a,b) convo_ranks(a, b), b = pmfs)
}



convo_lis1 <- function(lis, dat, vec)  {
  mclapply(seq_along(lis), function(a)
    mclapply(seq_along(a[[1]]), function(b)
      lapply(b, function(c)
        convo(dat[,lis[[a]][[b]][c]],dat[,lis[[a]][[b]][1]],vec))) %>% unlist) %>% setDT
}


#' derived pmf by subtraction
#'
#' Given two pmfs and a value vector, returns convolved pmf and value vector in matrix.
#' Values at associated with probs along (y) are subtracted from values associated with probs along (x)
#' throughout a discretized x,y event space.  Returns the set of prob:value pairs in a matrix
#' with cols prob, value ordered by value.
#'
#' @param x,y are sparse matrices with cols c(prob,value)
#' @return the set of prob:value pairs in a matrix with cols prob, value ordered by value.


convo_minus <- function(x,y)  {
  mat <- matrix(0,nrow=nrow(x)*nrow(y),ncol=2)
  mat[,1] <- mapply(function(a) mapply(function(b,c) b*c, b=a, c=y[,1]), a=x[,1])
  mat[,2] <- mapply(function(a) mapply(function(b,c) b-c, b=a, c=y[,2]), a=x[,2])
  mat[order(mat[,2]),]
}



#' convo_plus convolves two pmfs by adding their values.
#'
#' Given two pmfs and a value vector, returns convolved pmf and value vector in matrix.
#' Values at associated with probs along (x) are added to values associated with probs along (y)
#' throughout a discretized x,y event space.  Returns the set of prob:value pairs in a matrix
#' with cols prob, value ordered by value.
#'
#' @param x,y are sparse matrices with cols c(prob,value)
#' @return the set of prob:value pairs in a matrix with cols prob, value ordered by value.


convo_plus <- function(x,y)  {
  mat <- matrix(0,nrow=nrow(x)*nrow(y),ncol=2)
  mat[,1] <- mapply(function(a) mapply(function(b,c) b*c, b=a, c=y[,1]), a=x[,1])
  mat[,2] <- mapply(function(a) mapply(function(b,c) b+c, b=a, c=y[,2]), a=x[,2])
  mat[order(mat[,2]),]
}



#' convolve by rank
#'
#' Given a vector of ranks, pulls the pmfs for each rank and subtracts
#' the youngest from each via convolution.
#'
#' @param ranks is a vector of ranks (element of rank list)
#' @param pmfs is the matrix of charcoal pmfs
#' @param vals is the vector of values (years) associated with pmfs
#' @param dat is the charcoal variable data.table
#' @return pmfs of ranks convolved by subtracting youngest from each
convo_rank <- function(ranks,
                       pmfs = char_pmfs,
                       vals = as.numeric(rownames(pmfs)),
                       dat = charcoal)  {
  ranks <- sort(ranks)
  young <- pmfs[, colnames(pmfs) == dat[rank == ranks[1], site_id]]
  ar <- array(0,c(nrow(pmfs), length(ranks)))
  for (i in seq_along(ranks))  {
    vec <- pmfs[, colnames(pmfs) == dat[rank == ranks[i], site_id]]
    ar[,i] <- convo(young, vec, vals)
  }
  ar
}



#' convolve vector of ranks
#'
#' subtracts the pmf with first rank in `vec` from pmfs with ranks in `vec`
#' via convolution
#'
#' @param vec is a list of numeric sorted ranks
#' @param pmfs is a matrix of pmfs rows with cols of obs sorted by rank
#' @return a matrix of pmfs convolved by subtracting pmf of lowest rank from pmfs
#'   with ranks in `vec`
#' @import magrittr
#' @export

convo_ranks <- function(vec, pmfs = char_pmfs)  {
  index <- as.numeric(rownames(pmfs))
  mclapply(seq_along(vec), function(a)
    mclapply(seq_along(a[[1]]), function(b)
      lapply(b, function(c)
        convo(pmfs[,vec[[a]][[b]][c]],pmfs[,vec[[a]][[b]][1]], index))) %>% unlist) %>% setDT
}



#' Draw from exponential decay distribution
#'
#' Draw from exponential decay distribution $y = x^k$
#'
#' @param n is the max value returned by the decay function (integer)
#' @param k is the exponent $x^k$ (numeric)
#' @return \code{n} draws from the exponential decay distribution.
#' @export

draw_exp_decay <- function(n,k)  {
  xs <- 1:n
  ys <- (n - xs)^k
  unlist(mapply(function(x,y) rep(x,y), x = xs, y=ys))
}




#' Find the PMF of values vector.
#'
#' Given a numeric vector of values \code{vals}, and an index of values \code{index},
#' find the pmf of \code{vals}
#'
#' @param vals is a numeric vector of values
#' @param index is a numeric index of values
#' @return a numeric vector containing the pmf of \code{vals} along \code{index}
#' @export

fit_pmf <- function(vals,index)  {
  cdf <- to_cdf(vals)
  pmf <- to_pmf(cdf)
  mat <- matrix(c(pmf,vals),ncol=2)
  probcomb(mat,index)
}



#' get_cis returns confidence intervals of a vector
#'
#' Given a vector of probabilities (x), returns a numeric vector with
#' values at 2.5\%, median and 97.5\% of the distribution (x)
#'
#' @param x is a numeric vector of probabilities
#' @return a numeric vector with values at 2.5\%, median and 97.5\% of the distribution (x)
#' @export

get_cis <- function(x, lwr=.0250, med=.5000, upr=.9750)  {
  cis <- vector(length=3, mode='numeric')
  vec <- staTools::cdf(x)
  k <- 1
  while (vec$y[k]<lwr) k <- k+1
  cis[1] <- vec$x[k]
  while (vec$y[k]<med) k <- k+1
  cis[2] <- vec$x[k]
  while (vec$y[k]<=upr & k<length(vec$y)) k <- k+1
  cis[3] <- vec$x[k]
  names(cis) <- c('lwr','med','upr')
  cis
}


#' normalize
#'
#' normalize a vector of probabilities
#'
#' @param probs is a numeric vector of probabilities
#' @return normalized vector of `probs`
#' @export
normalize <- function(probs)  {
  norm <- array(0, length(probs))
  for (i in seq_along(probs))  {
    norm[i] <- probs[i] / sum(probs)
  }
  norm
}



#' pmf by rank
#'
#' wrapper for lapply, grabs pmfs with rank in `ranks` from main table
#'
#' @param ranks is a vector of ranks (element of rank list)
#' @param pmfs is the matrix of charcoal pmfs
#' @return the matrix of pmfs of rank `ranks`
#' @export

pmf_by_rank <- function(ranks, pmfs = char_pmfs)  {
  lapply(ranks, function(a,b) b[a,], b = pmfs)
}



#' Sums probability age value pairs into a single pmf
#'
#' Given a matrix of probability and age value pairs (x), and a vector of years (y),
#' returns a pmfs of length(y)
#'
#' @param x is a matrix of probability and age pairs
#' @param y is a numeric vector of years
#' @return a numeric vector of length(years) with value of sum probability age is less than
#'    year i in y and greater than year i-1
#'
#' @export


probcomb <- function(x,y)   {
  prob <- array(0,length(y))
  k <- 1
  n <- length(y)
  for (j in 1:nrow(x)){
    while (x[j,2] > y[k] & k < n)  k <- k + 1
    prob[k] <- prob[k] + x[j,1]
  }
  prob
}


#' Rank List
#'
#' Utility for returning sorted ranks of obs at sites above a minimum count.
#'
#' @param ftype is a character vector specifying facies type
#' @param dat is the data.table of charcoal site data
#' @param min_count is the minimum number of obs per site (1 is all sites)
#' @return a list of sorted ranks at sites with obs > `min_count`
#' @export

rank_list <- function(ftype, dat = charcoal, min_count = 3)  {
  dt <- dat[facies == ftype]
  dt_n <- dt[, .N]
  sites <- dt[, .N, keyby = .(family)]
  sites <- sites[N > min_count]
  lapply(unlist(sites[,1]), function(x) dt[family == x, rank])
}




#' sum_pmfs add pmfs together and returns a normalized sum
#'
#' Given a list of numeric pmfs (lis) and an integer representing the length
#' of pmfs in lis (len), returns a pmf representing the normalized sum of pmfs
#' in lis
#'
#' @param lis is a list of numeric pmfs
#' @param len is an integer representing the length of pmfs in lis
#' @return a summed and normalized pmf


sum_pmfs <- function(lis,len)  {
  mat <- matrix(unlist(lis), nrow = len)
  mat <- mapply(function(a,b,c) sum(b[a,]), a = seq_len(len), MoreArgs = list(b = mat))
  mat / sum(mat)
}




#' to_cdf converts numeric pmf vector to cdf of equal length
#'
#' Given a numeric pmf, returns a numeric vector of the cdf.
#'
#' @param vec is a numeric pmf
#' @return a numeric vector of the cdf
#'
#' @export

to_cdf <- function(vec)  {
  vec <- sort(vec)
  cdf <- array(0,length(vec))
  for (i in seq_along(vec)) {
    cdf[i] <- mean(vec<=vec[i])
  }
  cdf
}



#' to_exceed converts numeric pmf vectors to log10 exceedance distributions of equal length
#'
#' Given a numeric pmf, returns a numeric vector of the log10 exceedance distribution
#'
#' @param vec is a numeric pmf
#' @return a numeric vector of the log10 exceedance distribution
#' @seealso to_cdf

to_exceed <- function(vec)  {
  log10(1 - to_cdf(vec))
}


#' to_exceed converts numeric pmf vectors to log10 exceedance distributions of equal length
#'
#' Given a numeric pmf, returns a numeric vector of the log10 exceedance distribution
#'
#' @param vec is a numeric pmf
#' @return a numeric vector of the log10 exceedance distribution
#' @seealso to_cdf

to_lexc <- function(vec)  {
  log10(1 - to_cdf(vec))
}



#' CONVERT VECTOR TO pmf
#'
#' Given the cdf as a numeric vector, return the derived pmf as a numeric vector of equal length.
#'
#' @param cdf is a numeric vector (a cdf)
#' @return the pmf derived from \code{cdf} as a numeric vector
#' @export

to_pmf <- function(cdf)  {
  pmf <- array(0,length(cdf))
  pmf[1] <- cdf[1]
  for (i in 2:length(cdf))  {
    pmf[i] <- cdf[i] - cdf[i-1]
  }
  pmf
}


#' truncate at zero
#'
#' truncate a pmf at index values below zero
#'
#' @param pmf is a numeric vector of probabilites
#' @param index is a numeric vector values associates with `pmf`
#' @return `pmf` with probs at vals below zero set to zero
#' @export
trunc_0 <- function(pmf, index)  {
  trunc <- pmf
  for (i in seq_along(index))  {
    if (index[i] < 0)  trunc[i] <- 0
  }
  trunc
}


#' truncate & normalize
#'
#' truncates `pmf` at `index` values below zero and renormalizes
#'
#' @param pmf is a numeric vector of probabilities
#' @param index is a numeric vector of values associated with `pmf`
#' @return `pmf` truncated at `index` values below zero and renormalized
#' @export
trunc_norm <- function(pmf, index)  {
  normalize(trunc_0(pmf, index))
}


















#' Inherited age convolver
#'
#' By oversample site, draw a bootstrap sample with replacement (n) times,
#' convolve the sample distributions by subtracted the youngest from the set,
#' and return the 2.5\%, median, 97.5\% and observed distributions across all sites.
#'
#' @param dat is a data.table with rows of obs and cols of vars
#' @param ftype is a character vector matching facies class
#' @param n is an integer for number of bootstraps
#' @param t is a vector of years
#' @param probs is a matrix with row pmfs and cols of obs ordered by rank
#' @return matrix with row pmfs and cols for 2.5\%, median, 97.5\% and observed distributions
#'
#' @export
#' @import data.table
#' @import parallel
#' @importFrom magrittr %>%

boot_convo <- function(dat=sites, ftype, n, t=years, probs=npmf){
  # subset pmfs by facies type
  dt <- dat[facies == ftype]
  # number of obs
  dt_n <- dt[, .N]
  # number of obs by oss
  fc_n <- dt[, .N, keyby = .(family)]
  # rank list of ob rank by oss
  rl <- lapply(unlist(fc_n[,1]), function(x) dt[family == x, rank])
  #
  csl <- mclapply(seq_along(rl), function(a,b) convo_lis(b[[a]], probs, t), b=rl)
  csdt <- array(0,c(length(t),dt_n)) %>% as.data.table
  cism <- array(0,c(12,length(t)))
  csvec <- c(0,unlist(cumsum(fc_n[,2])))
  for (i in 1:length(csl))  {
    csdt[,(csvec[i]+1):csvec[i+1]] <- csl[[i]]
  }
  cism[4,] <- apply(csdt,1,function(x) sum(x)/ncol(csdt))
  cism[8,] <- cumsum(cism[4,])
  cism[12,] <- 1 - cism[8,]
  # id list of bootstrap samples of oss ranks by oss family, boot
  idl <- mclapply(rl, function(a) lapply(seq_len(n), function(x,y) sort(sample(y,rep=T)), y = a))

  cbl <- mclapply(seq_along(idl), function(a,b) convo_lis(b[[a]], probs, t), b=idl)
  dt <- array(0,c(length(t),length(cbl)*n)) %>% as.data.table
  for (i in 1:length(cbl))  {
    dt[,((i-1)*n+1):(i*n)] <- cbl[[i]]
  }
  cism[1:3,] <- apply(dt,1,get_cis)
  dt <- apply(dt,1,cumsum)
  bit <- apply(dt,2,get_cis)
  # dt <- apply(dt,1,function(a) 1-a)
  #  bat <- apply(dt,1,get_cis)
  csl
  #  cis <- mapply(function(a,b) get_cis(b[,.(a)]), a=seq_len(ncol(dt)), MoreArgs=list(b=dt))
  #  dt2 <- mapply(function(a,b) cumsum(b[,.(a)]/sum(b[,.(a)])), a=seq_len(ncol(dt)), MoreArgs=list(b=dt))
  #  cism[5:7,] <- mapply(function(a,b) get_cis(b[,.(a)]), a=seq_len(ncol(dt)), MoreArgs=list(b=dt))
  #  dt3 <- mapply(function(a,b) 1 - b[,.(a)], a=seq_len(ncol(dt2)), MoreArgs=list(b=dt2))
  #  cism[9:11,] <- mapply(function(a,b) get_cis(b[,.(a)]), a=seq_len(ncol(dt)), MoreArgs=list(b=dt))
  #  cis
}



#  DIFFER  #

differ <- function(vec,delta=0)	{

  #given a vector > length 2
  #return a variable of differences between elements

  for (i in 2:length(vec))	{
    delta[i] <- (vec[i] - vec[i-1])
  }
  return(delta)
}


# CUMULT #

cumult <- function(vec,dif,cum=0) {

  # given a vector > length 2
  # return a cumulative distribution curve

  for (i in 1:length(vec))	{
    cum[i] <- sum(vec[1:i]) / dif
  }
  return(cum)
}


#' Rack List
#'
#' Rack list with elements of equal length into a matrix using cbind
#'
#' @param lis is a list with elements of equal length
#' @return elements of `lis` column bound into a matrix
#' @export
rack <- function(lis)  {
  dt <- lis[[1]]
  for (i in 2:length(lis))  {
    dt <- cbind(dt, lis[[i]])
  }
  dt
}


#' CIs of Boot
#'
#' Return confidence intervals from a matrix of pmf bootstraps.
#'
#' @param boot is a matrix of pmfs
#' @param pmfs is the reference matrix of pmfs (charcoal pmfs)
#' @return vector of cis of `boot`
#' @import magrittr
#' @import parallel
#' @export
cis_of_boot <- function(boot, pmfs = char_pmfs)  {
  index <- as.numeric(rownames(pmfs)) %>% sort
  trunc <- mclapply(boot, function(a,b)  trunc_norm(a, b), b = index) %>% setDT
  dt <- trunc %>% t %>% as.data.table
  ob_pmf <- trunc %>% rowSums %>% normalize
  cis <- mclapply(dt, get_cis) %>% rack %>% rbind(ob_pmf)
  rownames(cis) <- c('lwr','med','upr','obs')
  colnames(cis) <- index %>% as.character
  cis
}

id_by_rank <- function(ranks, dat = charcoal)  {
  lapply(ranks, function(a) dat[rank == a, site_id]) %>% unlist
}






