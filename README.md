# ImportanceSamplingCode
Code to accompany 'Importance Sampling and Bayesian model comparison in ecology and evolution'

## Explantion of data required and format
(n is the number of individuals)

cint = a matrix [n, 2]; column 1 is either the last seen alive time (for interval censored individuals) or a zero (for right censored individuals). Column 2 is the recovered dead time (for interval censored individuals) or the last seen alive time (for right censored individuals)

censored = a vector of length n; this indicates the type of censoring for each individual. 1 = interval-censored, 2 = right-censored.

nind = number of individuals

sex = a vector of length n; this indicates the sex of each individual. 0 = female, 1 = male, NA = Unknown.

tD = a vector of length n; this is a dummy variable used by NIMBLE, all entries are NA.

zL = a vector of length n; this is the last seen alive time for all individuals.

zU = a vector of length n; this is either the recovered dead time (for interval censored individuals) or NA for right censored individuals.

