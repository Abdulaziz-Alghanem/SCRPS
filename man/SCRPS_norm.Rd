% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCRPS.R
\name{SCRPS_norm}
\alias{SCRPS_norm}
\title{SCRPS for the normal distribution}
\usage{
SCRPS_norm(y, mean = 0, sd = 1)
}
\arguments{
\item{y}{Vector of observations.}

\item{mean}{Vector of mean/location parameters.}

\item{sd}{Vector of standard deviation/scale parameters. Must be positive.}
}
\value{
Vector of scores.
}
\description{
Calculates the SCRPS for observations \code{y} and normal distributions with means = \code{mean} and standard deviatoins = \code{sd}.
}
\examples{
SCRPS_norm(0,mean = 0, sd = 1)
}
