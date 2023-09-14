###########################################################
##                                                      ###
##        Fit Gompertz model                            ###
##                                                      ###
###########################################################

## load libraries
library(nimble)
library(tidyverse)
library(mvtnorm)
library(boot)
library(lamW)
library(GGally)
library(coda)
library(mclust)
library(parallel)

rm(list=ls())

## set seed
set.seed(44)

## source necessary R functions
source("Dist_Gompertz.R")
source("Dist_GompertzNim.R")

## load data
load("CMR/BadgerData.RData")

## source additional R functions
source("ModelComparison_FUNCTIONS.R")

code <- nimbleCode({
  ## survival components for dead badgers
  for (i in 1:nind) {
    ## likelihood for interval-truncated gompertz
    censored[i] ~ dinterval(tD[i], cint[i, ])
    tD[i] ~ dgompzNim(a, b)
    
    ## sampling component
    pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), tM[i]) - y[i]) * log(1 - mean.p))
    dind[i] ~ dbern(pd[i])
  }
  
  ## priors
  a ~ dexp(1)
  b ~ dexp(1)
  mean.p ~ dunif(0,1)
})

## set up other components of model
consts <- list(nind = nind, tM = tM)
data <- list(y = y, cint = cint, 
             censored = censored, tD = tD, dind = dind)

## find overdispersed initial values
tinitFn <- function(cint, censored) {
  apply(cbind(cint, censored), 1, function(x) {
    if(x[3] == 2) {
      y <- x[2] + rexp(1, 1)
    } else {
      y <- runif(1, x[1], x[2])
    }
    y
  })
}
initFn <- function(cint, censored) {
  ## get ML estimates as initial values
  optFn <- function(pars, t) {
    if(any(pars < 0)) {
      return(NA)
    }
    sum(dGompertz(t, a = pars[1], b = pars[2], log = TRUE))
  }
  pars <- list(convergence = 1)
  k <- 0
  while(pars$convergence != 0 & k < 20) {
    ## sample missing values
    tD <- tinitFn(cint, censored)
    pars <- optim(rexp(2, 10), optFn, t = tD, control = list(fnscale = -1))
    k <- k + 1
  }
  if(k == 20) {
    stop("Can't sample initial values")
  }
  pars <- pars$par
  list(
    tD = tD,
    a = pars[1], 
    b = pars[2],
    mean.p = runif(1, 0, 1)
  )
}

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data, inits = initFn(cint, censored))

## compile the model
cmodel <- compileNimble(model)

## try with adaptive slice sampler
config <- configureMCMC(cmodel, monitors = c("a", "b", "mean.p"), thin = 1)
config$removeSamplers(c("a", "b"))
config$addSampler(target = c("a", "b"), type = 'AF_slice')

## check monitors and samplers
config$printMonitors()
config$printSamplers(c("a", "b", "mean.p"))

## build the model
built <- buildMCMC(config)
cbuilt <- compileNimble(built)

## run the model
system.time(run <- runMCMC(cbuilt, 
                             niter = 5000, 
                             nburnin = 190, 
                             nchains = 2, 
                             progressBar = TRUE, 
                             summary = TRUE, 
                             samplesAsCodaMCMC = TRUE, 
                             thin = 1))
run$summary

## plot mcmcm
plot(run$samples)
samples <- as.matrix(run$samples)

## pairs plot
samples <- as.matrix(samples)
samples <- samples[sample.int(nrow(samples), ceiling(nrow(samples) * 0.1)), ]
samples %>%
  as.data.frame() %>%
  ggpairs()

## fit range of finite mixture models
mod <- densityMclust(samples)

## summary of finite mixture models
summary(mod)
plot(mod, what = "BIC")

## take random samples from mixture
nimp <- 10000
nmix <- rbinom(1, size = nimp, prob = 0.95)
props <- sim(mod$modelName, mod$parameters, nmix)
props <- props[, -1, drop = FALSE]
colnames(props) <- c("a", "b", "mean.p")

## take random samples from prior (to create defense mixture)
dmp <- runif((nimp - nmix), 0, 1)
defense <- as.data.frame(matrix(rexp(2 * (nimp - nmix), 1), ncol = 2))
defense <- defense %>%
  mutate(mean.p = dmp)
colnames(defense) <- c("a", "b", "mean.p")

## check IS distribution against posterior samples
as.data.frame(props) %>%
  mutate(type = "IS") %>%
  rbind(as.data.frame(samples) %>%
          mutate(type = "Post")) %>%
  ggpairs(mapping = aes(colour = type, alpha = 0.5), upper = list(continuous = "density"), columns = 1:3)

## combine defense and importance samples
props <- rbind(props, defense)

## calculate log likelihood
loglike <- function(y, zL, zU, censored, tM, a, b, p) {
  ## loop over individuals
  fsurv1 <- pmap_dbl(list(y, zL, zU, censored, tM), function(y, zL, zU, censored, tM, a, b, p) {
    if(censored == 1){                             #interval
      ll <- y * log(p) + (zL - y) * log(1 - p)
      t <- (zL + 1):zU
      temp <- pGompertz(t, a, b) - pGompertz(t - 1, a, b)
      temp1 <- log(temp) + (t - 1 - zL) * log(1 - p)
      ll <- ll + log_sum_exp_marg(temp1, mn = FALSE)
    } else {                                     #right censored
      ll <- y * log(p) + (zL - y) * log(1 - p)
      if(zL < tM) {
        t <- (zL + 1):tM
        temp <- pGompertz(t, a, b) - pGompertz(t - 1, a, b)
        temp1 <- log(temp) + (t - 1 - zL) * log(1 - p)
        ## tail component (after last capture time)
        temp <- pGompertz(tM, a, b, lower.tail = FALSE)
        temp <- log(temp) + (tM - zL) * log(1 - p)
        temp1 <- c(temp1, temp)
        ll <- ll + log_sum_exp_marg(temp1, mn = FALSE)
      } else {
        temp <- pGompertz(tM, a, b, lower.tail = FALSE)
        ll <- ll + log(temp) + (tM - zL) * log(1 - p)
      }
    }
    ll
  }, a = a, b = b, p = p)
  sum(fsurv1)
}

## calculate log-likelihoods in parallel
logimpweight <- apply(props, 1, list)
logimpweight <- purrr::map(logimpweight, 1)
logimpweight <- mclapply(logimpweight,
                           function(pars, y, zL, zU, censored, tM) {
                             loglike(y, zL, zU, censored, tM, pars[1], pars[2], pars[3])
                           }, y = y, zL = zL, zU = zU, censored = censored, tM = tM, mc.cores = 24)
logimpweight <- reduce(logimpweight, base::c)

## add prior densities
logimpweight <- logimpweight + dexp(props[, 1], 1, log = TRUE) + 
  dexp(props[, 2], 1, log = TRUE) + dunif(props[, 3], 0, 1, log = TRUE)

## importance distributions
logimpweight <- logimpweight - 
  log(0.95 * dens(props, mod$modelName, mod$parameters, FALSE) + 0.05 * exp(dexp(props[, 1], 1, log = TRUE) +
                                                                              dexp(props[, 2], 1, log = TRUE) +
                                                                              dunif(props[, 3], 0, 1, log = TRUE)))

saveRDS(logimpweight, "CMR/outputs/logimpweight_g.rds")

## final checks
summary(props[is.finite(logimpweight), ])
summary(props)

## bootstrap the importance weights to create 95% intervals
imp.boot <- BootsPlot(logimpweight, 5000, TRUE)
