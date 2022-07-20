# library(INLA)
# library(SCRPS)
# library(doParallel)
# library(foreach)
#
# cs <- makeCluster(11)
# registerDoParallel(cs)




SCRPS_inla <- function(inla_result, n = 10000) {

  if (inla_result$.args$family == "gaussian") {

    return(SCRPS_inla_norm(inla_result = inla_result , n = n))

  } else if (inla_result$.args$family == "gamma") {


    return(SCRPS_inla_gamma(inla_result = inla_result , n = n))

  }







}



SCRPS_inla_norm <- function(inla_result, n = 10000) {

  i <- 0
  score = foreach::`%dopar%`(foreach::foreach(i=1:nrow(inla_result$summary.fitted.values)),  {

    mu_sample <- INLA::inla.rmarginal(n, inla_result$marginals.fitted.values[[i]])
    theta_sample <- INLA::inla.rmarginal(n,inla_result$marginals.hyperpar[[1]])
    y <- inla_result$.args$data$y[i]
    temp <- SCRPS::SCRPS_norm(y, mean = mu_sample, sd = 1/sqrt(theta_sample))
    score = mean(temp)
    return(score)
  })


  scores = unlist(score)
  return(mean(scores))


}



SCRPS_inla_gamma <- function(inla_result, n = 10000) {
  s <- rep(inla_result$.args$scale, length.out = nrow(inla_result$summary.fitted.values))
  i <- 0

  score = foreach::`%dopar%`(foreach::foreach(i = 1:nrow(inla_result$summary.fitted.values)), {

    mu_sample <- INLA::inla.rmarginal(n, inla_result$marginals.fitted.values[[i]])
    phi_sample <- INLA::inla.rmarginal(n,inla_result$marginals.hyperpar$`Precision parameter for the Gamma observations`)

    shape_sample <- phi_sample * s[i]
    scale_sample <- mu_sample / shape_sample
    y <- inla_result$.args$data$y[i]
    temp <- SCRPS::SCRPS_gamma(y, shape = shape_sample, scale = scale_sample)
    score = mean(temp)
    return(score)
  })


  scores = unlist(score)
  return(mean(scores))


}








