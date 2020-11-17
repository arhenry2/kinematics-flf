# New R Script to fit my groung truth parameters
## Copied fit_flf() from flf.R, but took it out of a function to get my fitted parameters from the ground truth
## Having them in a folder using the masterPath & res loading isn't working, so I'll have to break down the steps to see what's wrong
devtools::load_all('/Users/ashleyhenry/flf')
library(flf)
library(ggplot2)
library(plotly)
library(tidyverse)
library(readr)
library(rowr)
require(reshape2)
library(dplyr)
library(naniar)

# Makes gTruthPosVel with correct column names for flf_fit()
pos <- seq(0, 1500, 1) # pos values for ground truth parameters
# fitted vel values for ground truth parameters
## Ground truth parameters from RIL1_1/RIL1--RIL1_003_4--!!!
vel = flf(pos, 658.9043,	1.0581177,	0.014531372,	2.9281828) 
gTruthPosVel = cbind(pos, vel) # makes matrix of pos & fitted vel values for g truth parameters
gTruthPosVel = data.frame(gTruthPosVel) # forces the matrix to be a dataframe
# gTruthPosVel_1 = read.csv("~/Desktop/RILPop_tmp/RIL1_1/RIL1--RIL1_005_3--/rawData.csv")
# gTruthPosVel_1 = gTruthPosVel_1 %>%
#   rename(
#     pos = X10.055,
#     vel = X.0
#   )

fit_flf_1 <- function(inputList){
  library(stats4)
  # number we're using to define half width of the uniform distributions used for restarts (used in vfH, kH, nH)
  ## Note: used as st dev for rnorm() below
  per = 0.1
  pos <- gTruthPosVel_1$pos
  vel <- gTruthPosVel_1$vel
  # sort vf from highest to lowest values
  svel <- sort(vel, decreasing = TRUE) 
  # grab highest 100 vf values
  vfi <- mean(svel[1:100]) 
  print(c("mean of highest vf values", vfi))
  # find percent of vfi
  vfH <- vfi*per 
  v0i <- mean(vel)
  print(c("mean of vel values", v0i))
  x0i <- mean(pos)
  print(c("mean of pos values", x0i))
  x0H <- x0i*per
  ni <- 1.2
  nH <- ni*per
  ki <- .008
  kH <- ki*per
  # added 11.19.2019 to see if number of iterations through mle changes the estimated parameter values
  ## default = 100, had at 500 on 11.19.2019
  # maxitr = 100 
  maxitr = 100
  # added 11.19.2019 to see if difference between values at interations changes the estimated parameter values
  ## default = 1e-8, had at 1e-20 on 11.19.2019
  # reltol = 1e-20
  reltol = 1e-100
  # reN = 5
  reN = 10
  ll <- function(x0,vf,k,n){
    velp <- flf(pos,x0,vf,k,n)
    res <- vel-velp # calc's residuals
    # -sum(log(dnorm(res,0,1))) # old way before 3.2.2020
    # -mean(log(dnorm(res, 0, 0.1))) # new way on 11.10.2020, gives probability of your residuals falling into a N dist (assuming N dist.)
    -log(mean(dnorm(res,0,1))) #sd = 0.001
    #-sum(log(dnorm(res,0,0.1))) # Tried new one on 3.2.2020
    # -mean(log(dnorm(res,0,1)))
  }
  curMax = 0
  curMaxret = 0
  
  for (i in 1:reN){
    # Define new range for x0i (add noise for where computer will start looking for max likelihood estimate for parameter x0)
    x0ii = rnorm(1, mean = x0i, sd = x0H) 
    # Ditto as above, but for vf
    vfii = rnorm(1, mean = vfi, sd = vfH) 
    # Ditto as above, but for k
    kii = rnorm(1, mean = ki, sd = kH)
    # Ditto as above, but for n
    nii = rnorm(1, mean = ni, sd = nH) 
    start = list(x0 = 658.9043,	vf = 1.0581177,	k = 0.014531372,	n = 2.9281828) # type in answer: 658.90430000   1.05811770   0.01453137   2.92818280
    # print(start)
    ret <- mle(ll, start, method = "BFGS", control = list(maxit = maxitr, reltol = reltol))
    curV <- ll(coef(ret)[[1]], coef(ret)[[2]], coef(ret)[[3]], coef(ret)[[4]])
    if (curV > curMax){
      curMax = curV
      curMaxret = ret
    }
  
  # ret <- mle(ll, start = list(x0 = x0i,vf=vfi, k=ki, n=ni), method = "BFGS", control = list(maxit = maxitr, reltol = reltol))
  # ll(coef(ret)[[1]], coef(ret)[[2]], coef(ret)[[3]], coef(ret)[[4]])
  
  # mle(ll, start = list(x0 = x0i,vf=vfi, k=ki, n=ni)) # Original mle function w/out testing parameters for outliers
  # mle(ll, start = list(x0 = x0i,vf=vfi, k=ki, n=ni), method = "BFGS", control = list(maxit = maxitr, reltol = reltol)) # New mle function that tests for outliers
  
  # New mle function that tests for outliers
  curMaxret = mle(ll, start = list(x0 = x0i, vf=vfi, k=ki, n=ni), method = "BFGS", control = list(maxit = maxitr, reltol = reltol))

  list(coeffs = curMaxret, mle = curMax)
  coeffs <- list(coeffs = curMaxret, mle = curMax)
  }
  print(c("coeffs:", curMaxret))
  return(curMaxret)
  coeffs
}

fitTryData <- fit_flf_1(gTruthPosVel_1)

velFit = flf(pos, fitTryData@coef[[1]], fitTryData@coef[[2]], fitTryData@coef[[3]], fitTryData@coef[[4]])

### Add fitted vel from fit_flf() above
gTruthPosVelFit = cbind(pos, vel, velFit)
gTruthPosVelFit = data.frame(gTruthPosVelFit)
data_2 = cbind(pos, velFit)
data_2 = data.frame(data_2)

# Overlay plots for fitted vel (black) and ground truth vel (blue)
ggplot(gTruthPosVel_1, aes(x = pos)) +
  geom_point(aes(y = vel), color = "blue") +
  geom_point(data = data_2, aes(y = velFit)) +
  xlab("Position (px)") +
  ylab("Velocity (px/frame)") +
  ggtitle("Ground Truth & Incorrectly Fitted Curve") +
  theme_bw()



gTruth
# "gTruth" gives the following parameters:
#  x0           vf            k            n 
# 658.90430000   1.05811770   0.01453137   2.92818280
ret
# "ret" gives the following parameters:
#  x0           vf            k            n 
# 999.79974068   1.10537089   0.01387277   7.81947478 

# ret round infinity
#  x0           vf            k            n 
# 5.000061e+02 1.115285e+00 7.448094e-03 8.708169e-01 













  