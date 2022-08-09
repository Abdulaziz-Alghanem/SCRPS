

test_that("SCRPS for normal distribution", {

  n <- 100000
  k <- 50
  count <- 0
  for (i in 1:k) {
    y1 <- runif(1, min = -1000, max = 1000)
    m <- runif(1, min = -100, max = 100)
    s  <- runif(1, min = 1, max = 50)

    score1 <- SCRPS_norm(y1,mean = m, sd = s)
    sample1 <- rnorm(n,mean = m, sd = s)
    E1 <- mean(abs(sample1 - y1))
    sample1 <- rnorm(n,mean = m, sd = s)
    sample2 <- rnorm(n,mean = m, sd = s)
    E2 <- mean(abs(sample1 - sample2))
    score2 <- -E1/E2 - 0.5*log(E2)

    if (abs((score1 - score2)/(score1+10^-8)) < 0.05) {
      count <- count + 1
    }

  }

  accuracy = count/k

  expect_gt(accuracy, 0.95)











  } )




test_that("SCRPS for exponential distribution", {

  n <- 100000
  k <- 50
  count <- 0
  for (i in 1:k) {
    y1 <- runif(1, min = -1000, max = 1000)
    l1 <- runif(1, min = 1, max = 10)

    score1 <- SCRPS_exp(y1,rate=l1)
    sample1 <- rexp(n,rate=l1)
    E1 <- mean(abs(sample1 - y1))
    sample1 <- rexp(n,rate=l1)
    sample2 <- rexp(n,rate=l1)
    E2 <- mean(abs(sample1 - sample2))
    score2 <- -E1/E2 - 0.5*log(E2)

    if (abs((score1 - score2)/(score1+10^-8)) < 0.05) {
      count <- count + 1
    }

  }



  accuracy = count/k

  expect_gt(accuracy, 0.95)


} )





test_that("SCRPS for uniform distribution", {

  n <- 100000
  k <- 50
  count <- 0
  for (i in 1:k) {
    y1 <- runif(1, min = -1000, max = 1000)
    a <- runif(1, min = -100, max = 100)
    b <- runif(1, min = a + 1, max = 150)
    score1 <- SCRPS_unif(y1,min = a, max = b)
    sample1 <- runif(n,min = a, max = b)
    E1 <- mean(abs(sample1 - y1))
    sample1 <- runif(n,min = a, max = b)
    sample2 <- runif(n,min = a, max = b)
    E2 <- mean(abs(sample1 - sample2))
    score2 <- -E1/E2 - 0.5*log(E2)

    if (abs((score1 - score2)/(score1+10^-8)) < 0.05) {
      count <- count + 1
    }

  }

  accuracy = count/k

  expect_gt(accuracy, 0.95)


} )



test_that("SCRPS for exponential distribution with flexible location and scale", {

  n <- 100000
  k <- 50
  count <- 0
  for (i in 1:k) {
    y1 <- runif(1, min = -1000, max = 1000)
    l1 <- runif(1, min = 1, max = 10)
    a1 <- runif(1, min = -100, max = 100)
    b1 <- runif(1, min = 1, max = 50)

    score1 <- SCRPS_exp2(y1,rate=l1, location = a1, scale = b1)
    sample1 <- rexp(n,rate=l1)
    sample1 <- a1 + b1*sample1
    E1 <- mean(abs(sample1 - y1))
    sample1 <- rexp(n,rate=l1)
    sample1 <- a1 + b1*sample1
    sample2 <- rexp(n,rate=l1)
    sample2 <- a1 + b1*sample2
    E2 <- mean(abs(sample1 - sample2))
    score2 <- -E1/E2 - 0.5*log(E2)


    if (abs((score1 - score2)/(score1+10^-8)) < 0.05) {
      count <- count + 1
    }

  }




  accuracy = count/k

  expect_gt(accuracy, 0.95)


} )




test_that("SCRPS for logistic distribution", {

  n <- 100000
  k <- 50
  count <- 0
  for (i in 1:k) {
    y1 <- runif(1, min = -1000, max = 1000)
    m <- runif(1, min = -100, max = 100)
    s  <- runif(1, min = 1, max = 50)

    score1 <- SCRPS_logis(y1,location = m, scale = s)
    sample1 <- rlogis(n,location = m, scale = s)
    E1 <- mean(abs(sample1 - y1))
    sample1 <- rlogis(n,location = m, scale = s)
    sample2 <- rlogis(n,location = m, scale = s)
    E2 <- mean(abs(sample1 - sample2))
    score2 <- -E1/E2 - 0.5*log(E2)

    if (abs((score1 - score2)/(score1+10^-8)) < 0.05) {
      count <- count + 1
    }

  }


  accuracy = count/k

  expect_gt(accuracy, 0.95)


} )



test_that("SCRPS for laplace distribution", {

  library("ExtDist")
  n <- 100000
  k <- 50
  count <- 0
  for (i in 1:k) {
    y1 <- runif(1, min = -1000, max = 1000)
    m <- runif(1, min = -100, max = 100)
    s  <- runif(1, min = 1, max = 50)

    score1 <- SCRPS_lapl(y1,mu = m, b = s)
    sample1 <- rLaplace(n,mu = m, b = s)
    E1 <- mean(abs(sample1 - y1))
    sample1 <- rLaplace(n,mu = m, b = s)
    sample2 <- rLaplace(n,mu = m, b = s)
    E2 <- mean(abs(sample1 - sample2))
    score2 <- -E1/E2 - 0.5*log(E2)

    if (abs((score1 - score2)/(score1+10^-8)) < 0.05) {
      count <- count + 1
    }

  }





  accuracy = count/k

  expect_gt(accuracy, 0.95)


} )



test_that("SCRPS for two-piece normal distribution", {

  library("fanplot")

  n <- 100000
  k <- 50
  count <- 0
  for (i in 1:k) {
    y1 <- runif(1, min = -1000, max = 1000)
    m <- runif(1, min = -100, max = 100)
    s1  <- runif(1, min = 1, max = 50)
    s2  <- runif(1, min = 1, max = 50)

    score1 <- SCRPS_2pnorm(y1,mu = m, sd1 = s1, sd2 = s2)
    sample1 <- rsplitnorm(n, mode = m, sd1 = s1, sd2 = s2)
    E1 <- mean(abs(sample1 - y1))
    sample1 <- rsplitnorm(n, mode = m, sd1 = s1, sd2 = s2)
    sample2 <- rsplitnorm(n, mode = m, sd1 = s1, sd2 = s2)
    E2 <- mean(abs(sample1 - sample2))
    score2 <- -E1/E2 - 0.5*log(E2)

    if (abs((score1 - score2)/(score1+10^-8)) < 0.05) {
      count <- count + 1
    }

  }





  accuracy = count/k

  expect_gt(accuracy, 0.95)


} )








test_that("SCRPS for gamma distribution", {



  n <- 100000
  k <- 50
  count <- 0
  for (i in 1:k) {
    y1 <- runif(1, min = -1000, max = 1000)
    a <- runif(1, min = 1, max = 20)
    b <- runif(1, min = 1, max = 20)

    score1 <- SCRPS_gamma(y1,shape=a, scale = b)
    sample1 <- rgamma(n,shape= a, scale = b)
    E1 <- mean(abs(sample1 - y1))
    sample1 <- rgamma(n,shape= a, scale = b)
    sample2 <- rgamma(n,shape= a, scale = b)
    E2 <- mean(abs(sample1 - sample2))
    score2 <- -E1/E2 - 0.5*log(E2)

    if (abs((score1 - score2)/(score1+10^-8)) < 0.05) {
      count <- count + 1
    }

  }





  accuracy = count/k

  expect_gt(accuracy, 0.95)


} )


























test_that("SCRPS for Student's t distribution", {

  library("ggdist")

  n <- 100000
  k <- 50
  count <- 0
  for (i in 1:k) {
    y1 <- runif(1, min = -1000, max = 1000)
    df <- runif(1, min = 1, max = 100)
    mu <- runif(1, min = -100, max = 100)
    s <- runif(1, min = 1, max = 20)

    score1 <- SCRPS_t(y1,df = df, location = mu, scale = s)
    sample1 <- rstudent_t(n, df = df, mu = mu, sigma = s)
    E1 <- mean(abs(sample1 - y1))
    sample1 <- rstudent_t(n, df = df, mu = mu, sigma = s)
    sample2 <- rstudent_t(n, df = df, mu = mu, sigma = s)
    E2 <- mean(abs(sample1 - sample2))
    score2 <- -E1/E2 - 0.5*log(E2)

    if (abs((score1 - score2)/(score1+10^-8)) < 0.05) {
      count <- count + 1
    }

  }





  accuracy = count/k

  expect_gt(accuracy, 0.95)


} )













test_that("SCRPS for lognormal distribution", {


  n <- 100000
  k <- 50
  count <- 0
  for (i in 1:k) {
    y1 <- runif(1, min = -1000, max = 1000)
    mu <- runif(1, min = -7, max = 7)
    s <- runif(1, min = 0.2, max = 2.5)

    score1 <- SCRPS_lnorm(y1,meanlog = mu, sdlog = s)
    sample1 <- rlnorm(n, meanlog = mu, sdlog = s)
    E1 <- mean(abs(sample1 - y1))
    sample1 <- rlnorm(n, meanlog = mu, sdlog = s)
    sample2 <- rlnorm(n, meanlog = mu, sdlog = s)
    E2 <- mean(abs(sample1 - sample2))
    score2 <- -E1/E2 - 0.5*log(E2)

    if (abs((score1 - score2)/(score1+10^-8)) < 0.05) {
      count <- count + 1
    }

  }





  accuracy = count/k

  expect_gt(accuracy, 0.95)


} )






test_that("SCRPS for beta distribution", {


  n <- 100000
  k <- 50
  count <- 0
  for (i in 1:k) {
    y1 <- runif(1, min = -20, max = 20)
    a <- runif(1, min = 0.05, max = 20)
    b <- runif(1, min = 0.05, max = 20)

    score1 <- SCRPS_beta(y1,shape1 = a, shape2 = b)
    sample1 <- rbeta(n, shape1 = a, shape2 = b)
    E1 <- mean(abs(sample1 - y1))
    sample1 <- rbeta(n, shape1 = a, shape2 = b)
    sample2 <- rbeta(n, shape1 = a, shape2 = b)
    E2 <- mean(abs(sample1 - sample2))
    score2 <- -E1/E2 - 0.5*log(E2)

    if (abs((score1 - score2)/(score1+10^-8)) < 0.05) {
      count <- count + 1
    }

  }





  accuracy = count/k

  expect_gt(accuracy, 0.95)


} )















