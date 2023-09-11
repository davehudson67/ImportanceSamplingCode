##Functions used in Model Comparisons:

## calculates log of the marginal likelihood using log-sum-exp trick
log_sum_exp_marg <- function(x, ind = NULL, mn = TRUE) {
  
  if(!is.null(ind)) {
    x <- x[ind]
  }
  
  ## extract length
  n <- length(x)
  
  ## extract non-finite values
  ## (it's OK to remove these here because we
  ## assume they coorespond to a likelihood of zero
  ## and are thus felt through 1 / n above
  ## in the marginal likelihood calculation)
  x <- x[is.finite(x)]
  if(length(x) == 0) {
    return(-Inf)
  }
  
  ## extract maximum of logged values
  mX <- max(x)
  
  ## return answer
  out <- mX + log(sum(exp(x - mX)))
  if(mn) {
    out <- out - log(n)
  }
  out
}


## create and plot importance distribution against posterior samples
createImpDist <- function(samples, covscale = 1.2, draw = TRUE, grid_expanse = 0.1, ngrid = 50) {
  
  ## check samples are in correct format
  if(!is.mcmc(samples) & !is.mcmc.list(samples)) {
    stop("'samples' not in 'mcmc' or 'mcmc.list' format")
  }
  samples <- as.matrix(samples)
  
  ## calculate importance distribution
  impmean <- apply(samples, 2, mean)
  impcov <- cov(samples)
  diag(impcov) <- diag(impcov) * covscale
  
  if(draw) {
    
    ## extract number of parameters
    npars <- ncol(samples)
    
    ## create range of values to evaluate contours
    pgrid <- list(NULL)
    for(i in 1:npars) {
      agrid <- range(samples[, i])
      agrid[1] <- agrid[1] * ifelse(agrid[1] < 0, 1 + grid_expanse, 1 - grid_expanse)
      agrid[2] <- agrid[2] * ifelse(agrid[2] < 0, 1 - grid_expanse, 1 + grid_expanse)
      agrid <- seq(agrid[1], agrid[2], length.out = ngrid)
      pgrid[[i]] <- agrid
    }
    names(pgrid) <- colnames(samples)
    
    ## plot
    if(npars > 1) {
      ## create pairs
      inds <- combn(1:npars, 2)
      
      ## set up plot
      par(mfrow = c(ceiling(ncol(inds) / 3), min(ncol(inds), 3)))
      for(i in 1:ncol(inds)) {
        ## create grid
        grid <- expand.grid(pgrid[inds[, i]])
        
        ## calculate importance distribution across grid
        zgrid <- matrix(apply(grid, 1, dmvnorm, mean = impmean[inds[, i]], 
              sigma = impcov[inds[, i], inds[, i]]), ngrid, ngrid)
        
        ## plot
        plot(samples[, inds[1, i]], samples[, inds[2, i]], 
             xlab = colnames(samples)[inds[1, i]], ylab = colnames(samples)[inds[2, i]], pch = 20)
        contour(pgrid[[inds[1, i]]], pgrid[[inds[2, i]]], zgrid, add = TRUE, col = "red", lwd = 2)
      }
    } else {
      ## plot
      hist(samples[, 1], xlab = colnames(samples)[1], ylab = "Density", freq = FALSE)
      lines(pgrid[[1]], dnorm(pgrid[[1]], mean = impmean, sd = sqrt(impcov)), col = "red", lwd = 2)
    }
  }
  
  ## return values
  list(mean = impmean, cov = impcov)
}


## create and plot bootstrap samples of importance weights with 95% CIs
BootsPlot<-function(impWeights, r, trace = FALSE, nsteps = 10) {
  if(!trace) {
      boot.iw <- boot(data = impWeights, statistic = log_sum_exp_marg, R = r)
      ci <- boot.ci(boot.iw, type = "basic")
      hist(boot.iw$t)
      abline(v=ci$basic[4:5], col="red")
      br <- list(logmarg=boot.iw$t0, LCI=ci$basic[4],UCI=ci$basic[5])
  } else {
    
    ntot <- length(impWeights)
    if(ntot < 10) {
      stop("Number of importance samples < 10")
    }
    nsamps <- seq(10, ntot, length.out = nsteps)
    nsamps <- unique(round(nsamps))
    out <- NULL
    for (i in 1:length(nsamps)) {
      br <- impWeights[1:nsamps[i]]
      br <- boot(data = br, statistic = log_sum_exp_marg, R = r)
      br <- boot.ci(br, type = "basic")
      br <- data.frame(logmarg = br$t0, LCI = br$basic[4], UCI = br$basic[5])
      out[[i]] <- br
    }
    names(out) <- nsamps
    out <- bind_rows(out, .id = "nsamp") %>%
      mutate(nsamp = as.numeric(nsamp))
    
    print(ggplot(out, aes(x=nsamp)) +
      geom_point(aes(y = logmarg)) +
      geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5))
  }
  br
}

## create marginal likelihood plots
MargLike.plot <- function(boots) {
  m <- purrr::map(boots, as.data.frame)
  br <- bind_rows(m, .id = "model")
  br$model <- as.factor(br$model)
  
  ggplot(br, aes(x=model, y=logmarg, col=model)) +
    geom_point() +
    geom_errorbar(aes(ymin=LCI, ymax=UCI))
} 
  




