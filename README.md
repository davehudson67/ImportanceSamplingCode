# ImportanceSamplingCode
Code to accompany 'Importance Sampling and Bayesian model comparison in ecology and evolution'

## How to run the analysis

The code files here demonstrate fitting a Gompertz model to the Mongoose data and estimating the marginal likelihood. It then allows the comparison of this with 3 other different mortality models that have already been run.

Download the files from this repo and set the working directory to that folder.
Open and run the "FitGompertzMongoose.R" code which will then lead you through the process of completing marginal likelihood estimation via MCMC and Importance Sampling. The final stage of this code will save the output into the output folder.
Open and run the "InitialModelComparisonPostMongoose.R" code which will then compare the 4 mortality models and plot output.

"mong.RData" is the formatted Mongoose data.

## Explanation of data format
(n is the number of individuals)

cint = a matrix [n, 2]; column 1 is either the last seen alive time (for interval censored individuals) or a zero (for right censored individuals). Column 2 is the recovered dead time (for interval censored individuals) or the last seen alive time (for right censored individuals)

censored = a vector of length n; this indicates the type of censoring for each individual. 1 = interval-censored, 2 = right-censored.

nind = number of individuals

sex = a vector of length n; this indicates the sex of each individual. 0 = female, 1 = male, NA = Unknown.

tD = a vector of length n; this is a dummy variable used by NIMBLE, all entries are NA.

zL = a vector of length n; this is the last seen alive time for all individuals.

zU = a vector of length n; this is either the recovered dead time (for interval censored individuals) or NA for right censored individuals.

