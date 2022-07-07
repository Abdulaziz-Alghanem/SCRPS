library(INLA)
library(MASS)
library(SCRPS)
library(doParallel)
library(foreach)

cs <- makeCluster(11)
registerDoParallel(cs)

data(cement)

samplesize = 1000

x1 = rexp(samplesize)
x2 = runif(samplesize)
y = 1 + x1 + x2 + rnorm(samplesize)

df = data.frame(y=y,x1=x1,x2=x2)

# result <- inla(y ~ x1 + x2 + x3 + x4, data = cement, control.compute=list(return.marginals.predictor=TRUE))

result <- inla(y ~ x1 + x2, data = df, control.compute=list(return.marginals.predictor=TRUE))
summary(result)




MC_test <- function(inla_result, n = 10000) {

  # score <- inla_result$.args$data$y
  tot_number<-nrow(inla_result$summary.fitted.values)
  score = foreach (i=1:tot_number) %dopar% {
    mu_sample1 <- INLA::inla.rmarginal(n, inla_result$marginals.fitted.values[[i]])
    mu_sample2 <- INLA::inla.rmarginal(n, inla_result$marginals.fitted.values[[i]])
    #mu_sample3 <- inla.rmarginal(n, inla_result$marginals.fitted.values[[i]])
    theta_sample1 <- INLA::inla.rmarginal(n,inla_result$marginals.hyperpar[[1]])
    theta_sample2 <- INLA::inla.rmarginal(n,inla_result$marginals.hyperpar[[1]])
    #theta_sample3 <- inla.rmarginal(n,inla_result$marginals.hyperpar[[1]])
    #y_sample1 <- rnorm(n, mean = 1, sd = 10)
    #y_sample2 <- rnorm(n, mean = 1, sd = 10)
    y_sample1 <- rnorm(n, mean = mu_sample1, sd = 1/sqrt(theta_sample1))
    y_sample2 <- rnorm(n, mean = mu_sample2, sd = 1/sqrt(theta_sample2))
    #y_sample3 <- rnorm(n, mean = mu_sample3, sd = 1/sqrt(theta_sample3))
    y <- inla_result$.args$data$y[i]
    E1 <- mean(abs(y_sample1-y))
    E2 <- mean(abs(y_sample1-y_sample2))
    score <- -E1/E2 - 0.5*log(E2)
    return(score)
  }

  # print(score)
  scores = unlist(score)
  return(list(scores = scores, mean_scores=mean(scores)))










}

start = Sys.time()
MC_test(result)
end = Sys.time()
end-start











numeric_Int_test <- function(inla_result) {
  #

  tot_number <- nrow(inla_result$summary.fitted.values)

  scores = foreach(i = 1:tot_number) %dopar% {

    min_mu <- inla_result$marginals.fitted.values[[i]][1,"x"]
    max_mu <- inla_result$marginals.fitted.values[[i]][nrow(inla_result$marginals.fitted.values[[i]]),"x"]

    mu_dens <- function(x){
      dens <- sapply(x, function(z){
        if(z<min_mu){
          return(0)
        } else if (z>max_mu){
          return(0)
        } else{
          return(spline(x=inla_result$marginals.fitted.values[[i]][,"x"], y=inla_result$marginals.fitted.values[[i]][,"y"], xout=z)$y)

        }
      })
      return(dens)
    }


    min_theta <- inla_result$marginals.hyperpar[[1]][1,"x"]
    max_theta <- inla_result$marginals.hyperpar[[1]][nrow(inla_result$marginals.hyperpar[[1]]),"x"]

    theta_dens <- function(x){
      dens <- sapply(x, function(z){
        if(z<min_theta){
          return(0)
        } else if (z>max_theta){
          return(0)
        } else{
          return(spline(inla_result$marginals.hyperpar[[1]][,"x"], inla_result$marginals.hyperpar[[1]][,"y"], xout=z)$y)

        }
      })
      return(dens)
    }


    temp_func <- function(mu, theta) {

      return(SCRPS::SCRPS_norm(inla_result$.args$data$y[i], mean = mu, sd = 1/sqrt(theta))*mu_dens(mu)*theta_dens(theta))


    }





    #theta_grid <- seq(from = min_theta, to = max_theta, by = Step_size)
    theta_grid <- inla_result$marginals.hyperpar[[1]][,"x"]

    Integral1 = numeric(length=length(theta_grid))

    for(j in 1:length(theta_grid)){
      Integral1[j] <- integrate(f = temp_func, lower = min_mu, upper = max_mu, theta=theta_grid[j])[[1]]
    }


    temp_func2 <- function(x){
      dens <- sapply(x, function(z){
        if(z<theta_grid[1]){
          return(0)
        } else if (z>theta_grid[length(theta_grid)]){
          return(0)
        } else{
          return(spline(x=theta_grid, y=Integral1, xout=z)$y)

        }
      })
      return(dens)
    }



    Integral2 = integrate(f = temp_func2, lower = min_theta, upper = max_theta)[[1]]
    return(Integral2)

  }


  scores_vec <- unlist(scores)

  print(scores_vec)

  return(mean(scores_vec))






}

begin <- Sys.time()
numeric_Int_test(result)
end <- Sys.time()
end - begin





analytic_test <- function(inla_result, n = 10000) {


  score = foreach(i=1:nrow(inla_result$summary.fitted.values)) %dopar% {

    mu_sample <- INLA::inla.rmarginal(n, inla_result$marginals.fitted.values[[i]])
    theta_sample <- INLA::inla.rmarginal(n,inla_result$marginals.hyperpar[[1]])
    y <- inla_result$.args$data$y[i]
    #temp <- SCRPS_norm(y, mean = 1, sd = 10)
    temp <- SCRPS::SCRPS_norm(y, mean = mu_sample, sd = 1/sqrt(theta_sample))
    #temp <- SCRPS_norm(y, mean = mu_sample, sd = 1/sqrt(inla_result$summary.hyperpar$mean))
    score = mean(temp)
    return(score)
  }


  scores = unlist(score)
  print(scores)
  return(mean(scores))


}




start = Sys.time()
analytic_test(result)
end = Sys.time()
end-start



