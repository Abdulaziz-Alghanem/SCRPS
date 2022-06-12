

SCRPS_norm <- function(y, mean = 0, sd = 1) {
  z <- (mean-y) /sd
  score <- -sqrt(pi)*dnorm(z) - sqrt(pi)*z/2 * (2*pnorm(z) - 1) - 0.5*log(2*s/sqrt(pi))
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




