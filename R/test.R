score <- inla_result$.args$data$y



  mu_sample1 <- inla.rmarginal(n, inla_result$marginals.fitted.values[[i]])
  mu_sample2 <- inla.rmarginal(n, inla_result$marginals.fitted.values[[i]])
  theta_sample1 <- inla.rmarginal(n,inla_result$marginals.hyperpar[[1]])
  theta_sample2 <- inla.rmarginal(n,inla_result$marginals.hyperpar[[1]])
  y_sample1 <- rnorm(n, mean = 1, sd = 10)
  y_sample2 <- rnorm(n, mean = 1, sd = 10)
  # y_sample1 <- rnorm(n, mean = mu_sample1, sd = 1/sqrt(theta_sample1))
  # y_sample2 <- rnorm(n, mean = mu_sample2, sd = 1/sqrt(theta_sample2))
  y <- inla_result$.args$data$y[i]
  E1 <- cumsum(abs(y_sample1-y))/(1:10000)
  E2 <- cumsum(abs(y_sample1-y_sample2))/(1:10000)
  score[i] = -E1/E2 - 0.5*log(E2)
