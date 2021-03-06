% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{build_mats}
\alias{build_mats}
\title{build_mats condenses a sparse matrix down to a dense one.}
\usage{
build_mats(x, y)
}
\arguments{
\item{x}{is a numeric vector of probabilities (a pmf)}

\item{y}{is a numeric vector of values associated with x}
}
\value{
a dense matrix (x,y) excluding nonzero probabilities in x
}
\description{
Given two vectors of equal length, returns a matrix where values of the first vector are nonzero.
}
