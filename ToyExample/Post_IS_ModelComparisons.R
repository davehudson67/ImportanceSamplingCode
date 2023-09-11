## source additional R functions
source("ModelComparison_FUNCTIONS.R")

## load data
load("Sim_Data.RData")

## set up plot output file
pdf("outputs/IS_ModelComparisons.pdf")

###########################################################
##                                                      ###
##          Now conduct model comparisons               ###
##                                                      ###
###########################################################

## load IS samples
logimpweight_h0 <- readRDS("outputs/logimpweight_h0.rds")
logimpweight_h1 <- readRDS("outputs/logimpweight_h1.rds")

## generate log marginal likelihoods
logmarg_h0 <- log_sum_exp_marg(logimpweight_h0)
logmarg_h1 <- log_sum_exp_marg(logimpweight_h1)

## bootstrap samples
imp_boot_h0 <- BootsPlot(logimpweight_h0, 5000)
imp_boot_h1 <- BootsPlot(logimpweight_h1, 5000)

## add prior model weights
priorp <- 1/2
p_h0 <- logmarg_h0 + log(priorp)
p_h1 <- logmarg_h1 + log(priorp)

p <- c(p_h0, p_h1)
pd <- log_sum_exp_marg(p, mn = FALSE)

## normalise
p <- p - pd
p <- exp(p)
p

## plot marginal likelihoods
mods <- list(
  H0 = imp_boot_h0, 
  H1 = imp_boot_h1
)

MargLike.plot(mods)

## which models within log(20) of best
#logmarg <- map_dbl(mods, "logmarg")
#bestind <- which(logmarg == max(logmarg))
#logmargLCI <- mods[[bestind]]$LCI
#logmarg <- map_dbl(mods, "UCI")
#logmarg <- logmarg[map_lgl(logmarg, ~ . >= logmargLCI - log(20))]
#logmarg

## turn graphics device off
dev.off()
