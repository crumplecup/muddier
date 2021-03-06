% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mod.R
\name{slope_to_elev}
\alias{slope_to_elev}
\title{slope_to_elev}
\usage{
slope_to_elev(slope, dist, base = 0)
}
\arguments{
\item{slope}{is a numeric vector of slopes at creek nodes}

\item{dist}{is a numeric vector of distances to mouth at creek nodes}

\item{base}{is the starting elevation}
}
\value{
a vector of elevations
}
\description{
Converts slopes to elevations based on distance from a base.
}
