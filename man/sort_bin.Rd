% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score.R
\name{sort_bin}
\alias{sort_bin}
\title{Sort and bin x within a specified tolerance, using hist().}
\usage{
sort_bin(x, tol = 1e-05, nmax = NULL)
}
\arguments{
\item{x}{vector of covariates.}

\item{tol}{numeric tolerance, minimum distance between neighbouring points,
to avoid singularities.}

\item{nmax}{if specified, overrides tol as maximal number of unique points.}
}
\value{
list with three elements. x_sort is sorted and binned x,
w is a vector of weights corresponding to the frequency of each bin,
order is a vector specifying the ordering of x into the binned values
sort_x.
}
\description{
Sort and bin x within a specified tolerance, using hist().
}
