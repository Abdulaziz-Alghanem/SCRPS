#' Title
#'
#'
#' @param x b
#' @param ... b
#'
#' @return b
#' @export
#'
#' @examples
#' #to be completed
scrps <- function(x, ...) {

  UseMethod("scrps")


}












#' Title scrps for inla
#'
#' @param inla_result b
#' @param n b
#' @param test_data b
#' @param full_y b
#' @param test_idx b
#' @param only_test_data b
#' @param parallelize b
#' @param n.cores b
#'
#' @return b
#' @export
#'
#' @examples
#' #to be completed
scrps.inla <- function(inla_result,scoring_rule = "scrps", n = 10000,
                            test_data = NULL,
                            full_y = NULL,
                            test_idx = NULL,
                            only_test_data = TRUE,
                            parallelize = FALSE,
                            n.cores = parallel::detectCores()-1) {


  #work with the right y
  if(!is.null(test_data) || !is.null(full_y)){
    if(!is.null(full_y)){
      if(is.null(test_idx)){
        stop("You should provide the indexes for the test data in test_idx!")
      }
      y_complete <- full_y
      if (only_test_data) {
        test_idx_complete <- test_idx
      } else {
        test_idx_complete <- which(!is.na(y_complete))
      }

    } else{
      y_complete <- inla_result$.args$data$y
      if (only_test_data) {
        test_idx_complete <- which(is.na(y_complete))
        y_complete[is.na(y_complete)] <- test_data
      } else {
        y_complete[is.na(y_complete)] <- test_data
        test_idx_complete <- which(!is.na(y_complete))
      }

    }
  } else{
    y_complete <- inla_result$.args$data$y
    test_idx_complete <- which(!is.na(y_complete))
  }


  #work with the right s
  s_full <- inla_result$.args$scale
  if (is.null(s_full)) {
    s <- rep(1, length.out = nrow(inla_result$summary.fitted.values))
  } else if (length(s_full)==1){
    s <- rep(s_full, length.out = nrow(inla_result$summary.fitted.values))
  } else{
    s <- s_full
  }
  rm(s_full)



  if(parallelize){
    cs <-  parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cs)
    i <- 0
    # = foreach::`%dopar%`(foreach::foreach(i = 1:nrow(inla_result$summary.fitted.values)), {
    score = foreach::`%dopar%`(foreach::foreach(i = test_idx_complete), {

      temp <- inla_score_for_i(inla_result = inla_result, y = y_complete[i], s = s[i], i = i, n = n, scoring_rule = scoring_rule)
      score = mean(temp)
      return(score)
    })


    scores = unlist(score)
    parallel::stopCluster(cs)
    return(mean(scores))
  } else{
    #score = lapply(1:nrow(inla_result$summary.fitted.values), function(i){
    #score = lapply(which(!is.na(y_complete)), function(i){
    score = lapply(test_idx_complete, function(i){

      temp <- inla_score_for_i(inla_result = inla_result, y = y_complete[i], s = s[i], i = i, n = n, scoring_rule = scoring_rule)
      score = mean(temp)
      return(score)
    })


    scores = unlist(score)
    return(mean(scores))
    #return(scores)
  }



}











inla_score_for_i <- function(inla_result, y, s, i, n, scoring_rule) {

  #################################################################
  ####################      "Gaussian"      #########################
  #################################################################
  if (inla_result$.args$family == "gaussian") {
    mu_sample <- INLA::inla.rmarginal(n, inla_result$marginals.fitted.values[[i]])
    tau_sample <- INLA::inla.rmarginal(n,inla_result$marginals.hyperpar$`Precision for the Gaussian observations`)

    sd_sample <- 1/sqrt(tau_sample * s)


    if (scoring_rule == "scrps") {
      temp <- SCRPS::SCRPS_norm(y, mean = mu_sample, sd = sd_sample)

    } else if (scoring_rule == "crps") {

      temp <- -scoringRules::crps_norm(y, mean = mu_sample, sd = sd_sample)

    }
    return(temp)



    #################################################################
    ####################      "Gamma"      #########################
    #################################################################
  } else if (inla_result$.args$family == "gamma") {


    mu_sample <- INLA::inla.rmarginal(n, inla_result$marginals.fitted.values[[i]])
    phi_sample <- INLA::inla.rmarginal(n,inla_result$marginals.hyperpar$`Precision parameter for the Gamma observations`)

    shape_sample <- phi_sample * s
    scale_sample <- mu_sample / shape_sample


    if (scoring_rule == "scrps") {
      temp <- SCRPS::SCRPS_gamma(y, shape = shape_sample, scale = scale_sample)

    } else if (scoring_rule == "crps") {

      temp <- -scoringRules::crps_gamma(y, shape = shape_sample, scale = scale_sample)

    }



    return(temp)





    #################################################################
    ####################      "Log Normal"      #########################
    #################################################################
  } else if (inla_result$.args$family == "lognormal") {


    mu_sample <- INLA::inla.rmarginal(n, inla_result$marginals.fitted.values[[i]])
    tau_sample <- INLA::inla.rmarginal(n,inla_result$marginals.hyperpar$`Precision for the lognormal observations`)

    sd_sample <- 1/sqrt(tau_sample)

    if (scoring_rule == "scrps") {
      temp <- SCRPS::SCRPS_lnorm(y, meanlog = mu_sample, sdlog = sd_sample)

    } else if (scoring_rule == "crps") {

      temp <- -scoringRules::crps_lnorm(y, meanlog = mu_sample, sdlog = sd_sample)

    }


    return(temp)






  }



}







