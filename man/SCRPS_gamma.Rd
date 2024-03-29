% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SCRPS.R
\name{SCRPS_gamma}
\alias{SCRPS_gamma}
\title{SCRPS for Gamma distribution}
\usage{
SCRPS_gamma(y, shape = 1, rate = 1, scale = 1/rate)
}
\arguments{
\item{y}{Vector of observations.}

\item{shape}{Vector of shape parameters. Must be positive.}

\item{rate}{An alternative way to specify the scale.}

\item{scale}{Vector of scale parameters. Must be positive.}
}
\value{
Vector of scores.
}
\description{
Calculates the SCRPS for observations \code{y} and gamma distribution with shape parameter \code{shape} and scale parameter \code{scale} or alternatively rate parameter \code{rate} = 1/\code{scale}.
}
\examples{
#' SCRPS_gamma(7, 6, 2, 3)
}
