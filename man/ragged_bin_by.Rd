% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{ragged_bin_by}
\alias{ragged_bin_by}
\title{ragged_bin_by}
\usage{
ragged_bin_by(vec, by = vec, bins = 10)
}
\arguments{
\item{vec}{is a numeric vector of data}

\item{by}{is a range of values to bin by}

\item{bins}{is the number of bins to use}
}
\value{
a list length `bins` of index numbers of `vec` divided along `by`
}
\description{
divides range of `by` into `bins` intervals
returns list length `bins`
elements of list are index numbers of `vec`
where values are in range of interval of `by`
}
