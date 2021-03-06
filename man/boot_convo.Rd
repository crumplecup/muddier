% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{boot_convo}
\alias{boot_convo}
\title{Inherited age convolver}
\usage{
boot_convo(dat = sites, ftype, n, t = years, probs = npmf)
}
\arguments{
\item{dat}{is a data.table with rows of obs and cols of vars}

\item{ftype}{is a character vector matching facies class}

\item{n}{is an integer for number of bootstraps}

\item{t}{is a vector of years}

\item{probs}{is a matrix with row pmfs and cols of obs ordered by rank}
}
\value{
matrix with row pmfs and cols for 2.5\%, median, 97.5\% and observed distributions
}
\description{
By oversample site, draw a bootstrap sample with replacement (n) times,
convolve the sample distributions by subtracted the youngest from the set,
and return the 2.5\%, median, 97.5\% and observed distributions across all sites.
}
