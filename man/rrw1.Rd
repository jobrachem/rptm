% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/onion_spline.R
\name{rrw1}
\alias{rrw1}
\title{Generates samples from a random walk of order 1.}
\usage{
rrw1(n, p, sd = 1)
}
\arguments{
\item{n}{Positive integer, the number of samples.}

\item{p}{Positive integer, the number of parameters.}

\item{sd}{The standard deviation of the random walk.}
}
\value{
A numeric matrix of dimension \verb{(n, p)}.
}
\description{
Assumes the starting condition to be zero.
}
\examples{
samples <- rrw1(1, p = 10)
}
