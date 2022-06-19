

SCRPS_norm <- function(y, mean = 0, sd = 1) {
  z <- (mean-y) /sd
  score <- -sqrt(pi)*dnorm(z) - sqrt(pi)*z/2 * (2*pnorm(z) - 1) - 0.5*log(2*sd/sqrt(pi))

  return(score)
}



SCRPS_exp <- function(y, l = 1){
  if (y <= 0) {
    E1 <- 1/l - y
  } else {
    E1 <- (2*exp(-l*y))/l + y - 1/l
  }
  E2 <- 1/l
  score <- -E1/E2 - 0.5*log(E2)
  return(score)
}

SCRPS_unif <- function(y, min = 0, max = 1){
  if (y <= min) {
    E1 <- 0.5*(min + max) - y
  } else if (y >= max) {
    E1 <- y - 0.5*(min + max)
  } else {
    E1 <- (y*(y-min-max)+0.5*(min^2 + max^2))/(max-min)
  }

  E2 <- (max - min)/(3)
  score <- -E1/E2 -0.5*log(E2)
  return(score)
}



SCRPS_exp2 <- function(y, l = 1, a = 0, b = 1){
  if (y <= a) {
    E1 <- a + b/l - y
  } else {
    E1 <- (b/l)*(2*exp(-l*((y-a)/b)) - 1) + y - a
  }
  E2 <- b/l
  score <- -E1/E2 - 0.5*log(E2)
  return(score)
}


SCRPS_logis <- function(y, mu = 0, sd = 1){
  y <- (y-mu)/sd
  E1 <- 2*log(1 + exp(y)) -y

  E2 <- 2
  score <- -E1/E2 - 0.5*log(E2)
  score <- score - 0.5 * log(sd)
  return(score)
}


SCRPS_lapl <- function(y, mu = 0, b = 1){
  y <- (y-mu)/b
  if (y <= 0) {
    E1 <- exp(y) -y
  } else {
    E1 <- exp(-y) +y
  }
  E2 <- 3/2
  score <- -E1/E2 - 0.5*log(E2)
  score <- score - 0.5*log(b)
  return(score)
}



SCRPS_2pnorm <- function(y, mu = 0, s1 = 1, s2 = 1){
  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
  erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
  y <- (y-mu)
  if (y <= 0) {
    E1 <- ( (-1 + 2*exp((-y^2)/(2*s1^2)))*sqrt(2/pi)*s1^2 + sqrt(2/pi)*s2^2 + y*s1 - y*s2 + 2*s1*y*erf(y/(sqrt(2)*s1))  ) /(s1 + s2)
  } else {
    E1 <- (sqrt(2/pi)*s1^2 - sqrt(2/pi)*s2^2 + 2*exp(-y^2/(2*s2^2))*sqrt(2/pi)*s2^2 + s1*y + s2*y*erf(y/(sqrt(2)*s2)) -s2*y*erfc(y/sqrt(2)*s2)) / (s1 + s2)
  }
  E2 <-  (-2*(-2 + sqrt(2))*s1^2 + 4*(-1 + sqrt(2))*s1*s2 -2*(-2 + sqrt(2))*s2^2) / (sqrt(pi)*(s1 + s2))
  score <- -E1/E2 - 0.5*log(E2)
  return(score)
}






