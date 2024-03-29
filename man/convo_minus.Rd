% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{convo_minus}
\alias{convo_minus}
\title{derived pmf by subtraction}
\usage{
convo_minus(x, y)
}
\arguments{
\item{x, y}{are sparse matrices with cols c(prob,value)}
}
\value{
the set of prob:value pairs in a matrix with cols prob, value ordered by value.
}
\description{
Given two pmfs and a value vector, returns convolved pmf and value vector in matrix.
Values at associated with probs along (y) are subtracted from values associated with probs along (x)
throughout a discretized x,y event space.  Returns the set of prob:value pairs in a matrix
with cols prob, value ordered by value.
}
