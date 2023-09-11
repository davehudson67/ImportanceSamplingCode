## Importance sampling
library(tidyverse)
library(GGally)
library(mvtnorm)
library(boot)
library(GGally)
library(coda)
library(mclust)
library(parallel)
library(boot)

pdf("outputs/h1_plots.pdf")

## load samples
mcmc.out <- readRDS("samples/mcmc.out_h1.rds")

## pairs plot
samples <- as.matrix(mcmc.out$samples)
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
colnames(props) <- c("beta0", "beta1", "sigma")

## take random samples from prior (to create defense mixture)
defense <- matrix(rnorm(2 * (nimp - nmix), 0, 1), ncol = 2)
defense <- cbind(defense, runif(nimp - nmix, 0, 100))
colnames(defense) <- c("beta0", "beta1", "sigma")

## check IS distribution against posterior samples
as.data.frame(props) %>%
  mutate(type = "IS") %>%
  rbind(as.data.frame(samples) %>%
          mutate(type = "Post")) %>%
  ggpairs(mapping = aes(colour = type, alpha = 0.5), upper = list(continuous = "density"), columns = 1:3)

## combine defense and importance samples
props <- rbind(props, defense)

## generate importance weights
## log-likelihood function
log.like <- function(y, x, beta0, beta1, sigma) {
  mu <- beta0 + beta1 * x
  sum(dnorm(y, mu, sigma, log = TRUE))
}

## calculate log-likelihoods in parallel
logimpweight <- apply(props, 1, list)
logimpweight <- purrr::map(logimpweight, 1)
logimpweight <- mclapply(logimpweight,
                         function(pars, y, x) {
                           log.like(y, x, pars[1], pars[2], pars[3])
                         }, y = y, x = x, mc.cores = 24)
logimpweight <- reduce(logimpweight, c)

## priors
logimpweight <- logimpweight + dnorm(props[, 1], 0, 1, log = TRUE) +
                               dnorm(props[, 2], 0, 1, log = TRUE) + 
                               dunif(props[, 3], 0, 100, log = TRUE)

## importance distributions
logimpweight <- logimpweight - 
  log(0.95 * dens(props, mod$modelName, mod$parameters, FALSE) + 0.05 * exp(dnorm(props[, 1], 0, 1, log = TRUE) +
                                                                            dnorm(props[, 2], 0, 1, log = TRUE) + 
                                                                            dunif(props[, 3], 0, 100, log = TRUE)))
saveRDS(logimpweight, "outputs/logimpweight_h1.rds")

## final checks
summary(props[is.finite(logimpweight), ])
summary(props)

## calculate log-marginal likelihood
logmarg <- log_sum_exp_marg(logimpweight)

## bootstrap the importance weights and create 95% intervals
BootsPlot(logimpweight, 5000, trace = TRUE)

## turn graphics device off
dev.off()
