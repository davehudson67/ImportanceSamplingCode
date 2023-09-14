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

## Assuming all models have been fitted.

## load IS samples
logimpweight_e <- readRDS("CMR/outputs/logimpweight_e.rds")
logimpweight_g <- readRDS("CMR/outputs/logimpweight_g.rds")
logimpweight_gm <- readRDS("CMR/outputs/logimpweight_gm.rds")
logimpweight_s <- readRDS("CMR/outputs/logimpweight_s.rds")

## generate log marginal likelihoods
logmarg_e <- log_sum_exp_marg(logimpweight_e)
logmarg_g <- log_sum_exp_marg(logimpweight_g)
logmarg_gm <- log_sum_exp_marg(logimpweight_gm)
logmarg_s <- log_sum_exp_marg(logimpweight_s)

## bootstrap samples
imp_boot_e <- BootsPlot(logimpweight_e, 5000)
imp_boot_g <- BootsPlot(logimpweight_g, 5000)
imp_boot_gm <- BootsPlot(logimpweight_gm, 5000)
imp_boot_s <- BootsPlot(logimpweight_s, 5000)

## add prior model weights
pe <- logmarg_e + log(1/4)
pgm <- logmarg_gm + log(1/4)
pg <- logmarg_g + log(1/4)
ps <- logmarg_s + log(1/4)
p <- c(pe, pgm, pg, ps)
pd <- log_sum_exp_marg(p, mn = FALSE)

## normalise
p <- p - pd
p <- exp(p)
p

## plot marginal likelihoods
mods <- list(
  Exponential = imp_boot_e, 
  GompertzMakeham = imp_boot_gm, 
  Gompertz = imp_boot_g,
  Siler = imp_boot_s 
)
MargLike.plot(mods)

## which models within log(20) of best
logmarg <- map_dbl(mods, "logmarg")
bestind <- which(logmarg == max(logmarg))
logmargLCI <- mods[[bestind]]$LCI
logmarg <- map_dbl(mods, "UCI")
logmarg <- logmarg[map_lgl(logmarg, ~ . >= logmargLCI - log(20))]
logmarg


## Look at posterior model probabilities...

## group all logimpweights together
samples <- t(rbind(logimpweight_s, logimpweight_e, logimpweight_g, logimpweight_gm))

## set up bootsrap function
bootsIS <- function(samples, nboot, nmodel){
  lists <- list()
  priorp <- 1/nmodel
  
  for (i in 1:nboot) {
    samp_ints <- sample(1:10000, 10000, replace = TRUE)
    temp_post <- samples[samp_ints, ]
    logmarg <- apply(temp_post, 2, log_sum_exp_marg)
    logmarg <- logmarg + log(priorp)
    pd <- log_sum_exp_marg(logmarg, mn = FALSE)
    p <- logmarg - pd
    p <- exp(p)
    lists[[i]] <- p
    
  }
  t <- as.data.frame(do.call(rbind, lists))
  return(t)
}

## run function
out <- bootsIS(samples, 2000, 4)

## create CIs
IS_ModProbs <- out %>%
  apply(2, function(x) {
    quantile(x, probs = c(0.025, 0.5, 0.975))
  }) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("LCI", "Median", "UCI")) %>%
  mutate(model = colnames(out))

## adjust model names
IS_ModProbs$model <- c("Siler", "Exponential", "Gompertz", "Gompertz-Makeham")

## plot
ggplot(IS_ModProbs,aes(x = model, y = Median)) +
  geom_point() +
  geom_errorbar(ymax = IS_ModProbs$UCI, ymin = IS_ModProbs$LCI) +
  theme(axis.text.x = element_text(angle = 90))