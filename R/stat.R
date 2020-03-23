#' bin_cdfs
#'
#' returns a dataframe of change in cdf of each var by bin
#'
#' @param num is a numeric vector
#' @param den is a numeric vector
#' @param by is a numeric vector
#' @param rho is a number representing density
#' @param bins is the number of bins to use
#' @return a dataframe of change in cdf of each var by bin
#' @export

bin_cdfs <- function(num, den, by, rho = 1, bins = 10)  {
  num_bins <- ragged_bin_by(num , by, bins)
  den_bins <- ragged_bin_by(den , by, bins)
  rng <- num_bins[[2]]
  num_bins <- num_bins[[1]]
  den_bins <- den_bins[[1]]
  num_n <- 0
  den_n <- 0
  for(i in 1:bins) {
    num_n[i] <- length(num_bins[[i]])
    den_n[i] <- length(den_bins[[i]])
  }
  ndif <- (1 / rho) * (num_n / den_n)
  return(list(ndif, num_n, den_n, rng))
}



#' bin_weights
#'
#' returns a dataframe of change in cdf of each var by bin
#'
#' @param num is a numeric vector
#' @param den is a numeric vector
#' @param by is a numeric vector
#' @param bins is the number of bins to use
#' @return a dataframe of change in cdf of each var by bin
#' @export

bin_weights <- function(num, den, by, bins = 10)  {
  num_cdf <- to_cdf(num[order(by)] / sum(num))
  den_cdf <- to_cdf(den[order(by)] / sum(den))
  by_cdf <- to_cdf(by[order(by)] / sum(by))
  mat <- as.matrix(data.frame(num_cdf, den_cdf, by_cdf))
  vec_len <- length(by)
  bin_len <- ceiling(vec_len / bins)
  vec_ids <- 1:length(by)
  bin_list <- list()
  ranges <- vector(bins, mode = 'numeric')
  xs <- vector(bins, mode = 'numeric')
  ys <- vector(bins, mode = 'numeric')
  for (i in 1:bins)  {
    if (i == 1)  {
      bin_list[[i]] <- vec_ids[vec_ids <= bin_len]
      bin <- mat[bin_list[[i]], ]
      ranges[i] <- bin[nrow(bin), 3]
      xs[i] <- bin[nrow(bin),1] - bin[1,1]
      ys[i] <- bin[nrow(bin),2] - bin[1,2]
    }

    if (i > 1)  {
      bin_list[[i]] <- vec_ids[vec_ids > (i-1) * bin_len &
                                 vec_ids <= i * bin_len]
      bin <- mat[bin_list[[i]], ]
      ranges[i] <- bin[nrow(bin), 3]
      xs[i] <- bin[nrow(bin),1] - bin[1,1]
      ys[i] <- bin[nrow(bin),2] - bin[1,2]
    }
  }
  df <- data.frame(x = xs, y = ys, rng = ranges)
  return(df)
}



#' cor_creeks
#'
#' prints a correlation plot of two variables colored by creek
#'
#' @param var1 is a numeric vector
#' @param var2 is a numeric vector length(`var1`)
#' @param creeks is a character vector of creek names
#' @param labs is a character vector providing two plot labels in order `var1` `var2`
#' @param creek_labs is a character vector providing creek labels
#' @param leg_pos is a character vector denoting legend position
#' @param a is the level of alpha transparency (0,1)
#' @return a scatter plot of `var1` vs. `var2` with a linear fit line
#' @export

cor_creeks <- function(var1, var2, creeks, creek_labs,
                       labs = c('var1', 'var2'),
                       cols = length(creek_labs),
                       leg_pos = 'topleft', a = .33) {
  pal <- get_palette(cols, a)
  plot(var1, var2, pch = 20,
       col = rgb(red=.001, green=.001, blue=.001, alpha=.001),
       xlab = labs[1], ylab = labs[2],
       main = paste0('slope = ', round(lm(var1 ~ var2)$coefficients[2], 2),
                     '    R^2 = ', round(summary(lm(var1 ~ var2))$adj.r.squared, 2)))
  for (i in 1:length(creek_labs)) {
    points(var1[creeks %in% creek_labs[i]],
           var2[creeks %in% creek_labs[i]],
           pch = 20, col = pal[i])
  }
  abline(lm(var2 ~ var1), lwd = 2.5, col = 'goldenrod')
  legend(leg_pos, legend = creek_labs, fill = pal)
}


#' cor_plot
#'
#' prints a correlation plot of two variables
#'
#' @param var1 is a numeric vector
#' @param var2 is a numeric vector length(`var1`)
#' @param labs is a character vector providing two plot labels in order `var1` `var2`
#' @return a scatter plot of `var1` vs. `var2` with a linear fit line
#' @export

cor_plot <- function(var1, var2, labs = c('var1', 'var2')) {
  plot(var1, var2, pch = 20,
       col = rgb(red=.306, green=.569, blue=.831, alpha=.33),
       xlab = labs[1], ylab = labs[2],
       main = paste0('slope = ', round(lm(var1 ~ var2)$coefficients[2], 2),
                     '    R^2 = ', round(summary(lm(var1 ~ var2))$adj.r.squared, 2)))
  abline(lm(var2 ~ var1), lwd = 2.5, col = 'goldenrod')
}


#' discrete_bin
#'
#' divide a variable of interest `vals`
#' into intervals of length `bin_len`
#' return ids of `vals` as a list
#'
#' @param vals is a numeric vector, the variable of interest
#' @param bin_len is the interval to divide `vals` by
#' @return a list of ids for `vals` divided by `bin_len` into bins
#' @export

discrete_bin <- function(vals, bin_len) {
  n <- ceiling(max(vals))
  bins <- ceiling(n / bin_len)

  marks <- bin_len
  for (i in 2:(bins-1)) {
    marks[i] <- marks[i-1] + bin_len
  }
  marks <- c(marks, n)
  bin_vals <- list()
  for (i in 1:bins) {
    if (i == 1)  {
      bin_vals[[i]] <- vals[vals <= marks[i]]
    }
    if (i > 1) {
      bin_vals[[i]] <- vals[vals > marks[i-1] &
                              vals <= marks[i]]
    }
  }
  return(bin_vals)
}


#' fit_bins
#'
#' for each bin in `bins` select pmfs in `dat`
#' where expected means in bin equals `ids` in `dat`
#' compute interarrival times
#' as expected mean of convolved difference of successive arrivals
#' fit to exponential distribution
#' return interarrival times, exp fit, med and 95% cis
#'
#' @param dat is a numeric matrix of pmfs (pmf, samples)
#' @param ids is a numeric vector of expected means of samples in `dat`
#' @param bins is a list of numeric matrices of pmfs (pmf, samples)
#' @param dist is the distribution used by `fitdistrplus::fitdist`
#' @param index is a numeric vector of values associated with pmf in `dat`
#' @return a list of interarrival times, exp fits, med and ci values
#' @export

fit_bins <- function(dat, ids, bins,
                     dist = 'exp', index = as.numeric(rownames(char_pmfs))) {
  times <- list()
  fits <- list()
  coefs <- array(0, c(length(bins), 4))
  for (i in 1:length(bins))  {
    mat <- dat[ , ids %in% bins[[i]]]
    n <- ncol(mat)
    mat_mns <- array(0, c(n-1, 1))
    mat_ar <- array(0, c(n-1, nrow(mat)))

    for(j in 1:(n-1)) {
      mat_ar[j,] <- convo(mat[,j+1], mat[,j], index)  # convolved difference in age
      mat_mns[j] <- weighted.mean(sort(index), mat_ar[j,])
    }

    times[[i]] <- mat_mns
    fits[[i]] <- fitdistrplus::fitdist(as.numeric(mat_mns), dist)
    if(dist == 'exp') {
      coefs[i,1] <- summary(fits[[i]])$estimate
      coefs[i,2:4] <- fitdistrplus::bootdist(fits[[i]])$CI
    }
  }
  return(list(times, fits, coefs))
}



#' get_form
#'
#' from colnames of `df`, make into lm formula
#' first var is dependent, other vars independent
#'
#' @param df is a data.frame of data to model
#' @return a formula made from colnames of `df`
#' @export

get_form <- function(df) {
  labs <- colnames(df)
  n <- length(labs)
  form <- paste0(labs[1], ' ~ ')
  if (n > 2) {
    for (i in 2:(n-1)) {
      form <- paste0(form, labs[i], ' + ')
    }
  }
  form <- paste0(form, labs[n])
  form
}


#' get_forms
#'
#' from the colnames in `df`
#' first col is depvar, others indepvars
#' return character vector of all unique formulas
#' for feeding to lm() or glm()
#'
#' @param df is a data.frame object
#' @return character vector of formula combinations
#' @export
get_forms <- function(df) {
  nms <- colnames(df)
  ind <- nms[-1]
  dep <- paste0(nms[1], ' ~ ')
  for (i in seq_along(ind))  {
    k <- 1
    if (i == 1) {
      comb <- gtools::combinations(length(ind), i, ind)
      for (j in seq_along(comb)) {
        if (j == 1) {
          forms <- paste0(dep, comb[j])
          k <- k + 1
        }
        if (j > 1)  {
          forms[k] <- paste0(dep, comb[j])
          k <- k + 1
        }
      }
    }
    if (i > 1) {
      comb <- gtools::combinations(length(ind), i, ind)
      for (j in 1:nrow(comb)) {
        for (m in 1:ncol(comb)) {
          if (m == 1) {
            lab <- comb[j, m]
          }
          if (m > 1)  {
            lab <- paste(lab, comb[j, m], sep = ' + ')
          }
        }
        forms[k] <- paste0(dep, lab)
        k <- k + 1
      }
    }
  }
  return(forms)
}


#' get_palette
#'
#' get a palette of colors at chosen transparency
#' colors are coral, hardwood, gold, forest, sky, ocean
#' violet, rose, crimson
#' choose a palette with a vector of names
#' get a random palette length `x` if `x` is a number
#'
#' @param x is a vector of color names or number
#' @param a is the alpha transparency (0,1)
#' @return a palette of rgb values
#' @export

get_palette <- function(x, a = .33) {
  col_names <- c('coral', 'hardwood', 'gold', 'forest', 'leaf', 'sky', 'ocean',
            'violet', 'rose', 'crimson', 'white', 'slate', 'charcoal', 'black')
  cols <- c(rgb(.921, .251, .203, a, 'coral'),
            rgb(.321, .196, .129, a, 'hardwood'),
            rgb(.812, .675, .0, a, 'gold'),
            rgb(.0, .361, .024, a, 'forest'),
            rgb(.561, .82, .459, a, 'leaf'),
            rgb(.0, .753, .78, a, 'sky'),
            rgb(.02, .0, .612, a, 'ocean'),
            rgb(.608, .0, .89, a, 'violet'),
            rgb(.839, .369, .471, a, 'rose'),
            rgb(.788, .0, .0, a, 'crimson'),
            rgb(1, 1, 1, a, 'white'),
            rgb(.666, .666, .666, a, 'slate'),
            rgb(.333, .333, .333, a, 'charcoal'),
            rgb(0, 0, 0, a, 'black'))
  df <- data.frame(col_names = col_names, cols = cols,
                   stringsAsFactors = F)
  palette <- 0
  if (inherits(x, 'character')) {
    for (i in seq_along(x)) {
      palette[i] <- df$cols[df$col_names == x[i]]
    }
  }
  if (inherits(x, 'numeric')) {
    palette <- df$cols[sample.int(length(cols), x, replace = TRUE)]
  }
  return(palette)
}


#' gof_tab
#'
#' `fits` is a list of list (dist type(bins))
#' with elements of bins holding output from `fitdistrplus::fitdist`
#' fit using distribution dist type
#' `dists` specifies labels for distribution types for gof table
#' table of gof stats include bin no, dist type, aic, bic, ks, ad, cvm
#'
#' @param fits is a list of lists (dist type(bins))
#' @param dists is a character vector of dist type labels
#' @return a table of gof stats
#' @export

gof_tab <- function(fits, dists) {
  bin <- paste0('bin_' , sort(rep(1:length(fits[[1]]), length(fits))))
  dist_col <- rep(dists, length(fits))
  ar <- array(0, c(length(fits), length(fits[[1]]), 5))
  for (i in seq_along(fits)) {
    for (j in seq_along(fits[[i]])) {
      gof <- fitdistrplus::gofstat(fits[[i]][[j]])
      ar[i, j, 1] <- gof$aic
      ar[i, j, 2] <- gof$bic
      ar[i, j, 3] <- gof$ks
      ar[i, j, 4] <- gof$ad
      ar[i, j, 5] <- gof$cvm
    }
  }
  tab <- ar[ , 1, ]
  for (i in 2:length(fits[[1]])) {
    tab <- rbind(tab, ar[ , 2, ])
  }
  rownames(tab) <- bin
  colnames(tab) <- c('aic', 'bic', 'ks', 'ad', 'cvm')
  df <- data.table::data.table(bin = bin,
                               dist = dist_col,
                               aic = tab[ , 1],
                               bic = tab[ , 2],
                               ks = tab[ , 3],
                               ad = tab[ , 4],
                               cvm = tab[ , 5])
  df
}


#' logit_stat
#'
#' print logit model summary given data frame
#' first var is depedent, other vars independent
#'
#' @param df is a data.frame of data to model
#' @return a glm model summary fit to `df` using family binomial
#' @export

logit_stat <- function(df) {
  form <- get_form(df)
  mod <- glm(form, data = df, family = binomial)
  summary(mod)
}


#' mod_stat
#'
#' print lm model summary given data frame
#' first var is depedent, other vars independent
#'
#' @param df is a data.frame of data to model
#' @return a lm model summary fit to `df`
#' @export

mod_stat <- function(df) {
  form <- get_form(df)
  mod <- lm(form, df)
  summary(mod)
}


#' plot_weight
#'
#' prints a plot of the weighting function
#' produced using `bin_weights(num, dem, by, bin_no)`
#'
#' @param num is the numeric vector to use in the numerator of the weight
#' @param den is the numeric vector to use in the denomintor of the weight
#' @param by is the numeric vector to take the cdf of `num` and `dem` with respect to
#' @param bin_no is the number of bins to use in `bin_weights()`
#' @param labs is a character vector length 2 (xlab, ylab)
#' @return a plot of the weighting function using `weight_by()`
#' @seealso weight_by
#' @export

plot_weight <- function(num, den, by, bin_no = 10, labs) {
  wt <- weight_by(num, den, by, bin_no)
  plot(sort(by), wt, type = 'l', lwd = 2, col = 'slateblue',
       xlab = labs[1], ylab = labs[2])
  abline(h = 1, lwd = 1.5, lty = 2)
}


#' stat_tab
#'
#' in `df`, first var is dependent, other vars independent
#' gets all formula combinations
#' returns a table of model stats
#' adj.r.squared, AIC and BIC
#'
#' @param df is a data.frame object
#' @return table of lm statistics from `df`
#' @export

stat_tab <- function(df, type = 'lm') {
  forms <- get_forms(df)
  tab <- array(0, c(length(forms), 3))
  for (i in seq_along(forms)) {
    if (type == 'lm') {
      mod <- lm(forms[i], df)
      tab[i, 1] <- summary(mod)$adj.r.squared
      tab[i, 2] <- AIC(mod)
      tab[i, 3] <- BIC(mod)
    }
    if (type == 'binom') {
      mod <- glm(forms[i], data = df, family = 'binomial')
      tab[i, 1] <- rsq::rsq(mod)
      tab[i, 2] <- mod$aic
    }
  }
  return(data.frame(forms = forms,
                    R2 = tab[ , 1],
                    AIC = tab[ , 2],
                    BIC = tab[ , 3]))
}


#' weight_by
#'
#' generates a weighting function
#' produced using `bin_weights(num, dem, by, bin_no)`
#'
#' @param num is the numeric vector to use in the numerator of the weight
#' @param den is the numeric vector to use in the denomintor of the weight
#' @param by is the numeric vector to take the cdf of `num` and `dem` with respect to
#' @param bin_no is the number of bins to use in `bin_weights()`
#' @return a numeric vector of weighting function
#' @seealso bin_weights
#' @export

weight_by <- function(num, den, by, bin_no = 10) {
  wt_df <- bin_weights(num, den, by, bin_no)
  wt_val <- wt_df$x / wt_df$y
  by_cdf <- to_cdf(by[order(by)] / sum(by))
  wt <- vector(length(by_cdf), mode = 'numeric')
  for (i in seq_along(wt_df$rng)) {
    if (i == 1) {
      wt[by_cdf <= wt_df$rng[i]] <- wt_val[i]
    }
    if (i > 1) {
      wt[by_cdf > wt_df$rng[i-1] & by_cdf <= wt_df$rng[i]] <- wt_val[i]
    }
  }
  return(wt)
}












