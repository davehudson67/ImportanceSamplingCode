## Bayesian model comparison techniques blog
## load libraries
library(tidyverse)
library(nimble)
library(bayestestR)

## source additional functions
source("ModelComparison_FUNCTIONS.R")

## set up linear regression data
# generate independent variable values
x <- seq(-5, 5, by = 0.1)

# generate true parameter values
beta_0_true <- 1.2
beta_1_true <- 0.6
sigma_true <- 1

# generate dependent variable values
y <- beta_0_true + beta_1_true * x + rnorm(length(x), mean = 0, sd = sigma_true)
plot(x, y)

save.image("Sim_Data.RData")
