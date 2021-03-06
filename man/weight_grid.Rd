% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mod.R
\name{weight_grid}
\alias{weight_grid}
\title{weight_grid}
\usage{
weight_grid(pos, xy_grid, rec, test, by, scl = 1)
}
\arguments{
\item{pos}{is the coordinate pair of data in the record}

\item{xy_grid}{is the topology array to replace}

\item{rec}{is the record of test values}

\item{test}{is a string naming the test type c('ks', 'kp)}

\item{by}{is a string naming the test varaible c('vol', 'lvl')}

\item{scl}{is a scaling parameter, a float}
}
\value{
the estimated test value at `pos`
}
\description{
Returns topology of coordinate on grid given record of test values.
Types for `test`:
Kolmogorov-Smirnov statistic = 'ks'
Kuiper statistic = 'kp'
`xy_grid` determines the size of the outgoing topology array.
Types for `by`:
Unit cross-sectional area in m2 = 'vol'
Average cross-sectional depth per node = 'lvl'
`scl` ranges [0,1] and deflates the effects of neighbors on each other.
}
