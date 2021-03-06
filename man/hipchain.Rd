% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spat.R
\name{hipchain}
\alias{hipchain}
\title{hipchain}
\usage{
hipchain(pt, dist, inc = 10, chan = flowline, us = T, filter = 4000)
}
\arguments{
\item{pt}{is a coordinate pair}

\item{dist}{is a numeric hipchain distance (meters along channel)}

\item{inc}{is a numeric step increment}

\item{chan}{is a matrix of coords delineating channel from downstream to upstream}

\item{us}{is a boolean indicating up or downstream travel}

\item{filter}{is a numeric buffer distance defining off-channel pts}
}
\value{
coordinate position `dist` along `chan` from `pt`, and distance remaining in list
}
\description{
Given a matrix representation of the channel `chan` with cols (x,y) and rows running
upstream to downstream, pulls all pts within `inc` of `pt` along `chan`.
If `us` is true, then the last row of the subset matrix is `inc` units upstream
If `us` is false, then the first row of the subset is `inc` units downstream
The function recursively calls itself using the new point position and
distance remaining in place of `pt` and `dist` until having fully traveled
the original `dist` by increments of `inc`.
}
