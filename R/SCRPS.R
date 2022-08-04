

#'@name SCRPS_norm
#'
#'@title SCRPS for the normal distribution
#'
#'@description
#'Calculates the SCRPS for observations \code{y} and normal distributions with means = \code{mean} and standard deviatoins = \code{sd}.
#'
#' @param y Vector of observations.
#' @param mean Vector of mean/location parameters.
#' @param sd Vector of standard deviation/scale parameters. Must be positive.
#'
#' @return Vector of scores.
#' @export
#'
#' @examples
#' SCRPS_norm(0,mean = 0, sd = 1)
SCRPS_norm <- function(y, mean = 0, sd = 1) {
  z <- (mean-y) /sd
  score <- -sqrt(pi)*stats::dnorm(z) - sqrt(pi)*z/2 * (2*stats::pnorm(z) - 1) - 0.5*log(2*sd/sqrt(pi))

  return(score)
}


#'@name SCRPS_exp
#'
#'@title SCRPS for the exponential distribution
#'
#'@description
#'Calculates the SCRPS for observations \code{y} and exponential distributions with rates = \code{rate}.
#'
#' @param y Vector of observations.
#' @param rate Vector of rate parameters. Must be positive.
#'
#' @return Vector of scores.
#' @export
#'
#' @examples
#' SCRPS_exp(2, rate = 3)
SCRPS_exp <- function(y, rate = 1){
  ifelse(y<=0,E1 <- 1/rate - y, E1 <- (2*exp(-rate*y))/rate + y - 1/rate )
  E2 <- 1/rate
  score <- -E1/E2 - 0.5*log(E2)
  return(score)
}



#'@name SCRPS_unif
#'
#'@title SCRPS for the uniform distribution
#'
#'@description
#'Calculates the SCRPS for observations \code{y} and uniform distributions with lower limits = \code{min} and upper limits = \code{max}.
#'
#' @param y Vector of observations.
#' @param min Vector of lower limits.
#' @param max Vector of upper limits.
#'
#' @return Vector of scores.
#' @export
#'
#' @examples
#' SCRPS_unif(3, min = 10, max = 17)
SCRPS_unif <- function(y, min = 0, max = 1){
  ifelse(y<=min,E1 <- 0.5*(min + max) - y, ifelse(y >= max, E1 <- y - 0.5*(min + max),E1 <- (y*(y-min-max)+0.5*(min^2 + max^2))/(max-min) ) )

  E2 <- (max - min)/(3)
  score <- -E1/E2 -0.5*log(E2)
  return(score)
}



#' SCRPS for exponential distribution with flexible location and scale
#'
#'@description
#'Calculates the SCRPS for observations \code{y} and exponential distributions with rate parameters = \code{rate}, location parameters = \code{location} and scale parameters = \code{scale}.
#'
#' @param y Vector of observations.
#' @param rate Vector of rate parameters. Must be positive.
#' @param location Vector of location parameters.
#' @param scale Vector of scale parameters. Must be positive.
#'
#' @return Vector of scores.
#' @export
#'
#' @examples
#' SCRPS_exp2(10, rate = 2, location = 5, scale = 3)
SCRPS_exp2 <- function(y, rate = 1, location = 0, scale = 1){
  ifelse(y <= location, E1 <- location + scale/rate - y,  E1 <- (scale/rate)*(2*exp(-rate*((y-location)/scale)) - 1) + y - location )
  E2 <- scale/rate
  score <- -E1/E2 - 0.5*log(E2)
  return(score)
}


#' SCRPS for logistic distribution
#'
#' @description
#' Calculates the SCRPS for observations \code{y} and logistic distributions with location parameters = \code{location} and scale parameters = \code{scale}.
#'
#' @param y Vector of observations.
#' @param location Vector of location parameters.
#' @param scale Vector of scale parameters. Must be positive.
#'
#' @return Vector of scores.
#' @export
#'
#' @examples
#' SCRPS_logis(3,location = -5, scale = 2)
SCRPS_logis <- function(y, location = 0, scale = 1){
  y <- (y-location)/scale
  E1 <- 2*log(1 + exp(y)) -y

  E2 <- 2
  score <- -E1/E2 - 0.5*log(E2)
  score <- score - 0.5 * log(scale)
  return(score)
}


#' SCRPS for the laplace distribution
#'
#' @description
#' Calculates the SCRPS for observations \code{y} and laplace distributions with location parameters = \code{mu} and scale parameters = \code{b}.
#'
#' @param y Vector of observations.
#' @param mu Vector of location parameters.
#' @param b Vector of scale parameters. Must be positive.
#'
#' @return Vector of scores.
#' @export
#'
#' @examples
#' SCRPS_lapl(2, mu = 3, b = 2)
SCRPS_lapl <- function(y, mu = 0, b = 1){
  y <- (y-mu)/b
  ifelse(y <= 0 ,  E1 <- exp(y) -y, E1 <- exp(-y) +y)
  E2 <- 3/2
  score <- -E1/E2 - 0.5*log(E2)
  score <- score - 0.5*log(b)
  return(score)
}



#' SCRPS for two-piece normal distribution
#'
#' @description
#'  Calculates the SCRPS for observations \code{y} and two-piece normal distributions with location parameters = \code{mu} and left-hand-side scale parameters = \code{sd1} and right-hand-side scale parameters = \code{sd2}.
#'
#' @param y Vector of observations.
#' @param mu Vector of location parameters.
#' @param sd1 Vector of left-hand-side scale parameters. Must be positive.
#' @param sd2 Vector of right-hand-side scale parameters. Must be positive.
#'
#' @return Vector of scores.
#' @export
#'
#' @examples
#' SCRPS_2pnorm(7, 6, 2, 3)
SCRPS_2pnorm <- function(y, mu = 0, sd1 = 1, sd2 = 1){
  erf <- function(x) 2 * stats::pnorm(x * sqrt(2)) - 1
  erfc <- function(x) 2 * stats::pnorm(x * sqrt(2), lower = FALSE)
  y <- (y-mu)
  ifelse(y <= 0,
         E1 <- ( (-1 + 2*exp((-y^2)/(2*sd1^2)))*sqrt(2/pi)*sd1^2 + sqrt(2/pi)*sd2^2 + y*sd1 - y*sd2 + 2*sd1*y*erf(y/(sqrt(2)*sd1))  ) /(sd1 + sd2),
         E1 <- (sqrt(2/pi)*sd1^2 - sqrt(2/pi)*sd2^2 + 2*exp(-y^2/(2*sd2^2))*sqrt(2/pi)*sd2^2 + sd1*y + sd2*y*erf(y/(sqrt(2)*sd2)) -sd2*y*erfc(y/sqrt(2)*sd2)) / (sd1 + sd2))
  E2 <-  (-2*(-2 + sqrt(2))*sd1^2 + 4*(-1 + sqrt(2))*sd1*sd2 -2*(-2 + sqrt(2))*sd2^2) / (sqrt(pi)*(sd1 + sd2))
  score <- -E1/E2 - 0.5*log(E2)
  return(score)
}







#' SCRPS for Gamma distribution
#'
#' @description
#'  Calculates the SCRPS for observations \code{y} and gamma distribution with shape parameter \code{shape} and scale parameter \code{scale} or alternatively rate parameter \code{rate} = 1/\code{scale}.
#'
#' @param y Vector of observations.
#' @param shape Vector of shape parameters. Must be positive.
#' @param rate An alternative way to specify the scale.
#' @param scale Vector of scale parameters. Must be positive.
#'
#' @return Vector of scores.
#' @export
#'
#'
#' @examples
#' #' SCRPS_gamma(7, 6, 2, 3)

SCRPS_gamma <- function(y, shape = 1, rate = 1, scale = 1/rate){
  E1 <- 2*y*stats::pgamma(y,shape = shape, scale = scale) - y - 2*shape*scale*stats::pgamma(y, shape = shape+1, scale = scale) + shape*scale
  E2 <- 2*scale/beta(1/2, shape)
  score <- -E1/E2 - 0.5*log(E2)
  return(score)
}










#' SCRPS for Student's t distribution
#'
#'@description
#'  Calculates the SCRPS for observations \code{y} and Student's t distribution with degrees of freedom parameter \code{df} and optionally location parameter \code{location} and scale parameter \code{scale}.
#'
#' @param y Vector of observations.
#' @param df Vector of degrees of freedom parameters. Must be greater than 1.
#' @param location Vector of location parameters.
#' @param scale Vector of scale parameters. Must be positive.
#'
#' @return Vector of scores.
#' @export
#'
#' @examples
#' SCRPS_t(1,2,3,4)
#'
SCRPS_t <- function(y, df, location = 0, scale = 1) {
  y <- (y - location)/scale

  E1 <- y*(2*stats::pt(y,df = df) - 1) + 2*stats::dt(y, df =df)*(df + y^2)/(df - 1)
  E2 <- (4*sqrt(df))*(beta(1/2, df- 1/2))/(df-1)/(beta(1/2, df/2))^2
  score <- -E1/E2 - 0.5*log(E2)
  score <- score - 0.5 * log(scale)
  return(score)
}










#' SCRPS for log-normal distribution
#'
#'@description
#'  Calculates the SCRPS for observations \code{y} and lognormal distribution with locationlog parameter \code{meanlog} and scalelog parameter \code{sdlog}.
#'
#' @param y Vector of observations.
#' @param meanlog Vector of locationlog parameters.
#' @param sdlog Vector of scalelog parameters. Must be positive.
#'
#' @return Vector of scores.
#' @export
#'
#' @examples
#' SCRPS_lnorm(1,2,3)
SCRPS_lnorm <- function(y, meanlog = 0, sdlog = 1){
  ifelse(y <= 0, temp <- 0,  temp <- stats::pnorm((log(y) - meanlog - sdlog^2)/sdlog))
  E1 <- y*(2*stats::plnorm(y,meanlog = meanlog, sdlog = sdlog) - 1) -2*exp(meanlog + (sdlog^2)/2)*(temp-0.5)
  E2 <- 4*exp(meanlog + (sdlog^2)/2)*(stats::pnorm(sdlog/sqrt(2))-0.5)
  score <- -E1/E2 - 0.5*log(E2)
  return(score)
}









#'
#'
#'
#' #' SCRPS for log-laplace distribution
#' #'
#' #'@description
#' #'  Calculates the SCRPS for observations \code{y} and loglaplace distribution with locationlog parameter \code{locationlog} and scalelog parameter \code{scalelog}.
#' #'
#' #' @param y Vector of observations.
#' #' @param locationlog Vector of locationlog parameters.
#' #' @param scalelog Vector of scalelog parameters. Must be between 0 and 1
#' #'
#' #' @return Vector of scores.
#' #' @export
#' #'
#' #' @examples
#' #' SCRPS_llapl(1,2,0.7)
#' SCRPS_llapl <- function(y, locationlog = 0, scalelog = 0.5){
#'
#'   FF <- function(x, locationlog, scalelog ) {
#'     ifelse(x <= 0,0,0.5*(1 + sign(log(x)- locationlog)*(1-exp(-abs(log(x)- locationlog)/scalelog))))
#'
#'   }
#'
#'   A <- function(x, locationlog, scalelog) {
#'
#'     ifelse(x < exp(locationlog), (1-(2*FF(x,locationlog,scalelog))^(1+scalelog))/(1+scalelog) , (1- (2*(1- FF(x,locationlog,scalelog)))^(1-scalelog))/(scalelog - 1) )
#'
#'
#'
#'   }
#'
#'
#'   E1 <- y*(2*FF(y,locationlog,scalelog) - 1) + exp(locationlog)*(A(y,  locationlog, scalelog))
#'   E2 <-  0.25*exp(locationlog)*((scalelog/(scalelog^2-4)))
#'   score <- -E1/E2 - 0.5*log(E2)
#'   return(score)
#' }
#'


# \eqn(x+y) <-> $x+y$
# \deqn(x+y) <-> $$x+y$$
