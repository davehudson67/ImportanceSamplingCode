#Create custom Gompertz distribution:
## probability density function
dgompzNim <- nimbleFunction(
  run = function(x = double(0), a = double(0),
                 b = double(0),
                 log = integer(0, default = 0)) {
    returnType(double(0))
    if(a < 0 | b < 0 ) {
      return(NaN)
    }
    logS <- (a / b) * (1 - exp(b * x))
    logH <- log(a) + (b * x)
    logProb <- logH + logS
    if(log) return(logProb)
    else return(exp(logProb))
  })

## function to produce random samples
rgompzNim <- nimbleFunction(
  run = function(n = integer(0), a = double(0),
                 b = double(0)) {
    returnType(double(0))
    if(a < 0 | b < 0 ) {
      return(NaN)
    }
    if(n != 1) print("rgompzNim only allows n = 1; using n = 1.")
    u <- runif(1, 0, 1)
    rs <- (1 / b) * log(1 - (b / a) * log(1 - u))
    return(rs)
  })

## cumulative distribution function (and survivor function)
pgompzNim <- nimbleFunction(
  run = function(q = double(0), a = double(0),
                 b = double(0),
                 lower.tail = integer(0, default = 1), 
                 log.p = integer(0, default = 0)) {
    returnType(double(0))
    if(a < 0 | b < 0 ) {
      return(NaN)
    }
    logS <- (a / b) * (1 - exp(b * q))
    if(!lower.tail) { 
      if(log.p) return(logS)
      else return(exp(logS))
    } else {
      p <- 1 - exp(logS)
      if(!log.p) return(p)
      else return(log(p))
    }
  })

## quantile function
qgompzNim <- nimbleFunction(
  run = function(p = double(0), a = double(0),
                 b = double(0),
                 lower.tail = integer(0, default = 1), 
                 log.p = integer(0, default = 0)) {
    returnType(double(0))
    if(a < 0 | b < 0 ) {
      return(NaN)
    }
    if(log.p) p <- exp(p)
    if(!lower.tail) p <- 1 - p
    return((1 / b) * log(1 - (b / a) * log(1 - p)))
  })

## register distributions with NIMBLE
registerDistributions(list(
  dgompzNim = list(
    BUGSdist = "dgompzNim(a, b)",
    Rdist = "dgompzNim(a, b)",
    pqAvail = TRUE, 
    range = c(0, Inf)
  )
))
