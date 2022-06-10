

SCRPS_norm <- function(y, mean = 0, sd = 1) {
  z <- (mean-y) /sd
  score <- -sqrt(pi)*dnorm(z) - sqrt(pi)*z/2 * (2*pnorm(z) - 1) - 0.5*log(2*s/sqrt(pi))
  return(score)
}


SCRPS_exp <- function(y, l = 1){
   score <- -2*exp(-l*y) - l*y + 1 - 0.5*log(1/l)
   return(score)
  }




