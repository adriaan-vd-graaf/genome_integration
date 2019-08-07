# Title     : This is copied from:
#           Consistent Estimation in Mendelian Randomization with Some Invalid Instruments Using a Weighted Median Estimator
#           Jack Bowden, 1 George Davey Smith, 1 Philip C. Haycock, 1 and Stephen Burgess
#           PMID: 27061298
# Objective : Does a weighted median analysis.
# Created by: adriaan
# Created on: 8/6/19

weighted.median <- function(betaIV.in, weights.in) {
    betaIV.order = betaIV.in[order(betaIV.in)]
    weights.order = weights.in[order(betaIV.in)]
    weights.sum = cumsum(weights.order)-0.5*weights.order
    weights.sum = weights.sum/sum(weights.order)
    below = max(which(weights.sum<0.5))
    weighted.est = betaIV.order[below] + (betaIV.order[below+1]-betaIV.order[below])* (0.5-weights.sum[below])/(weights.sum[below+1]-weights.sum[below])
    return(weighted.est)
}

weighted.median.boot <- function(betaXG.in, betaYG.in, sebetaXG.in, sebetaYG.in, weights.in){
    med = NULL
    for(i in 1:1000){
        betaXG.boot = rnorm(length(betaXG.in), mean=betaXG.in, sd=sebetaXG.in)
        betaYG.boot = rnorm(length(betaYG.in), mean=betaYG.in, sd=sebetaYG.in)
        betaIV.boot = betaYG.boot/betaXG.boot
        med[i] = weighted.median(betaIV.boot, weights.in)
    }
    return(sd(med))
}

args = commandArgs(trailingOnly=TRUE)
# cat(args, "\n")

table = read.table(args[1], sep="\t", header=TRUE)
# print(table)

betaXG = table[,1]
sebetaXG = table[,2]
betaYG = table[,3]
sebetaYG = table[,4]

betaIV = betaYG/betaXG

# ratio estimates
weights = (sebetaYG/betaXG)^-2
# inverse-variance weights
betaIVW = sum(betaYG*betaXG*sebetaYG^-2)/sum(betaXG^2*sebetaYG^-2)
# penalized weights
betaWM = weighted.median(betaIV, weights)
cat(betaWM, "\n")
# weighted median estimate
sebetaWM = weighted.median.boot(betaXG, betaYG, sebetaXG, sebetaYG, weights)
cat(sebetaWM, "\n")
# IVW estimate
# standard error
# penalized weighted median estimate