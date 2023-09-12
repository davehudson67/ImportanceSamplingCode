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
library(survival)
rm(list = ls())

## source necessary R distributions
source("Dist_Gompertz.R")
source("Dist_GompertzNim.R")

## source additional R functions
source("ModelComparison_FUNCTIONS.R")

## load data
load("mong.RData")

## set seed according to model
set.seed(2)

###########################################################
##                                                      ###
##        Now fit Gompertz model                        ###
##                                                      ###
###########################################################

code <- nimbleCode({
  ## survival components for dead badgers
  for (i in 1:nind) {
  ## likelihood for interval-truncated gompertz makeham
    censored[i] ~ dinterval(tD[i], cint[i, ])
    tD[i] ~ dgompzNim(a, b)
    
  }
  
  ## priors
  a ~ dexp(1)
  b ~ dexp(1)
})

## set up other components of model
consts <- list(nind = nind)
data <- list(cint = cint, censored = censored, tD = tD)

tinitFn <- function(cint, censored) {
  apply(cbind(cint, censored), 1, function(x) {
    if(x[3] == 2) {
      y <- x[2] + runif(1, 0, 1)
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
    pars <- optim(rexp(2, 100), optFn, t = tD, control = list(fnscale = -1))
    k <- k + 1
  }
  if(k == 20) {
    stop("Can't sample initial values")
  }
  pars <- pars$par
  list(
    tD = tD,
    a = pars[1], 
    b = pars[2]
  )
}

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data, inits = initFn(cint, censored))

## compile the model
cmodel <- compileNimble(model)

## set monitor
config <- configureMCMC(cmodel, monitors = c("a", "b"), thin = 1)
config$removeSamplers(c("a", "b"))
config$addSampler(target = c("a", "b"), type = 'AF_slice')

## check monitors and samplers
config$printMonitors()
config$printSamplers(c("a", "b"))

## build the model
built <- buildMCMC(config)
cbuilt <- compileNimble(built)

## run the model
system.time(run <- runMCMC(cbuilt, 
    niter = 20000, 
    nburnin = 5000, 
    nchains = 2, 
    progressBar = TRUE, 
    summary = TRUE, 
    samplesAsCodaMCMC = TRUE, 
    thin = 1))

## plot mcmc
plot(run$samples)
samples <- run$samples

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
props <- props[, -1]
colnames(props) <- c("a", "b")

## take random samples from prior (to create defense mixture)
defense <- matrix(rexp(2 * (nimp - nmix), 1), ncol = 2)
colnames(defense) <- c("a", "b")

## check IS distribution against posterior samples
as.data.frame(props) %>%
  mutate(type = "IS") %>%
  rbind(as.data.frame(samples) %>%
          mutate(type = "Post")) %>%
  ggpairs(mapping = aes(colour = type, alpha = 0.5), upper = list(continuous = "density"), columns = 1:2) +
  theme(text = element_text(size = 16))

## combine defense and importance samples
props <- rbind(props, defense)

## generate importance weights
## log-likelihood function
stopifnot(length(table(censored)) == 2)
log.like <- function(zL, zU, censored, a, b) {
  ## calculate log-likelihoods
  llI <- log(pGompertz(zU[censored == 1], a, b) - pGompertz(zL[censored == 1], a, b))
  llR <- pGompertz(zL[censored == 2], a, b, log = TRUE, lower.tail = FALSE)
  sum(llI) + sum(llR)
}

## calculate log-likelihoods in parallel
logimpweight <- apply(props, 1, list)
logimpweight <- purrr::map(logimpweight, 1)
logimpweight <- mclapply(logimpweight,
    function(pars, zL, zU, censored) {
       log.like(zL, zU, censored, pars[1], pars[2])
    }, zL = zL, zU = zU, censored = censored, mc.cores = 24)
logimpweight <- reduce(logimpweight, c)

## priors
logimpweight <- logimpweight + dexp(props[, 1], 1, log = TRUE) + 
  dexp(props[, 2], 1, log = TRUE)

## importance distributions
logimpweight <- logimpweight - 
  log(0.95 * dens(data = props, modelName = mod$modelName, parameters = mod$parameters, logarithm = FALSE) + 
        0.05 * exp(dexp(props[, 1], 1, log = TRUE) + 
     dexp(props[, 2], 1, log = TRUE)))
saveRDS(logimpweight, "outputs/logimpweight_g.rds")

## final checks
summary(props[is.finite(logimpweight), ])
summary(props)

## calculate log-marginal likelihood
logmarg <- log_sum_exp_marg(logimpweight)

## bootstrap the importance weights and create 95% intervals
BootsPlot(logimpweight, 5000, trace = TRUE)
