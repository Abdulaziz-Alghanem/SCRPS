

#' SCRPS for INLA objects
#'
#' @param inla_result inla result object
#' @param n bla
#' @param test_data ble
#' @param full_y bli
#' @param test_idx blo
#' @param parallelize blu
#' @param n.cores bla
#'
#' @return SCRPS
#' @export
#'
#' @examples
#'
#' print("To be completed")

SCRPS_inla <- function(inla_result, n = 10000,
                       test_data = NULL,
                       full_y = NULL,
                       test_idx = NULL,
                       parallelize = FALSE,
                       n.cores = parallel::detectCores()-1) {

  if (inla_result$.args$family == "gaussian") {

    return(SCRPS_inla_norm(inla_result = inla_result , n = n,
                           test_data,
                           full_y,
                           test_idx,
                           parallelize = parallelize,
                           n.cores = n.cores))

  } else if (inla_result$.args$family == "gamma") {


    return(SCRPS_inla_gamma(inla_result = inla_result , n = n,
                            test_data,
                            full_y,
                            test_idx,
                            parallelize = parallelize,
                            n.cores = n.cores))

  }







}



SCRPS_inla_norm <- function(inla_result, n = 10000,
                            test_data,
                            full_y,
                            test_idx,
                            parallelize,
                            n.cores) {



  if(!is.null(test_data) || !is.null(full_y)){
    if(!is.null(full_y)){
      if(is.null(test_idx)){
        stop("You should provide the indexes for the test data in test_idx!")
      }
      y_complete <- full_y
    } else{
      y_complete <- inla_result$.args$data$y
      y_complete[is.na(y_complete)] <- test_data
    }
  }

  if(parallelize){
    cs <- parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cs)
    i <- 0
    score = foreach::`%dopar%`(foreach::foreach(i=1:nrow(inla_result$summary.fitted.values)),  {

      mu_sample <- INLA::inla.rmarginal(n, inla_result$marginals.fitted.values[[i]])
      theta_sample <- INLA::inla.rmarginal(n,inla_result$marginals.hyperpar[[1]])
      y <- y_complete[i]


      temp <- SCRPS::SCRPS_norm(y, mean = mu_sample, sd = 1/sqrt(theta_sample))
      score = mean(temp)
      return(score)
    })


    scores = unlist(score)
    parallel::stopCluster(cs)
    return(mean(scores))
  } else{
    score = lapply(1:nrow(inla_result$summary.fitted.values), function(i){

      mu_sample <- INLA::inla.rmarginal(n, inla_result$marginals.fitted.values[[i]])
      theta_sample <- INLA::inla.rmarginal(n,inla_result$marginals.hyperpar[[1]])
      y <- y_complete[i]


      temp <- SCRPS::SCRPS_norm(y, mean = mu_sample, sd = 1/sqrt(theta_sample))
      score = mean(temp)
      return(score)
    })


    scores = unlist(score)
    return(mean(scores))
  }



}



SCRPS_inla_gamma <- function(inla_result, n = 10000,
                             test_data,
                             full_y,
                             test_idx,
                             parallelize,
                             n.cores) {

  if(!is.null(test_data) || !is.null(full_y)){
    if(!is.null(full_y)){
      if(is.null(test_idx)){
        stop("You should provide the indexes for the test data in test_idx!")
      }
      y_complete <- full_y
    } else{
      y_complete <- inla_result$.args$data$y
      y_complete[is.na(y_complete)] <- test_data
    }
  }

  s_full <- inla_result$.args$scale
  if(length(s_full)==1){
    s <- rep(s_full, length.out = nrow(inla_result$summary.fitted.values))
  } else{
    s <- s_full
    rm(s_full)
  }



  if(parallelize){
    cs <-  parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cs)
    i <- 0
    score = foreach::`%dopar%`(foreach::foreach(i = 1:nrow(inla_result$summary.fitted.values)), {

      mu_sample <- INLA::inla.rmarginal(n, inla_result$marginals.fitted.values[[i]])
      phi_sample <- INLA::inla.rmarginal(n,inla_result$marginals.hyperpar$`Precision parameter for the Gamma observations`)

      shape_sample <- phi_sample * s[i]
      scale_sample <- mu_sample / shape_sample
      y <- y_complete[i]

      temp <- SCRPS::SCRPS_gamma(y, shape = shape_sample, scale = scale_sample)
      score = mean(temp)
      return(score)
    })


    scores = unlist(score)
    parallel::stopCluster(cs)
    return(mean(scores))
  } else{
    score = lapply(1:nrow(inla_result$summary.fitted.values), function(i){

      mu_sample <- INLA::inla.rmarginal(n, inla_result$marginals.fitted.values[[i]])
      phi_sample <- INLA::inla.rmarginal(n,inla_result$marginals.hyperpar$`Precision parameter for the Gamma observations`)

      shape_sample <- phi_sample * s[i]
      scale_sample <- mu_sample / shape_sample
      y <- y_complete[i]

      temp <- SCRPS::SCRPS_gamma(y, shape = shape_sample, scale = scale_sample)
      score = mean(temp)
      return(score)
    })


    scores = unlist(score)
    return(mean(scores))
  }



}








