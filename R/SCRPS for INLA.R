library(INLA)
library(MASS)
source("SCRPS.R")
data(cement)

print(cement)


result <- inla(y ~ x1 + x2 + x3 + x4, data = cement, control.compute=list(return.marginals.predictor=TRUE))
summary(result)

print(result$summary.fitted.values)
print(result$summary.linear.predictor)
print(result$summary.fixed)
print(result$summary.hyperpar)
SCRPS_norm(1,mean = 0, sd = 1)
SCRPS_LR <- function(inla_result) {
















}



