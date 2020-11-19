# install.packages("neldermead")
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
library(neldermead)

masterPath <- "/Users/ashleyhenry/Desktop/RILPop_full"

res <- processMasterFolder(masterPath)

#############################################################################################################################
# Below is code taken from flf.R, where I overlaid plots from nls() and Nelder-Mead via fminsearch() to find the best fitting
## Saved a few plots on Desktop to show proof of Nelder-Mead fitting better (will put into powerpoint for later)
### This was done on 11.18.2020
#############################################################################################################################

# Load in rawData to try different examples:
posVel = read.csv("~/Desktop/RILPop_full/RIL150_1/RILpop--RIL150_1--RIL150_001_3--/rawData.csv", header = FALSE)

####################################################################################
# Nelder-Mead Method
####################################################################################
# Tried with fminsearch(), which uses neldermead & is more closely related to matlab's fitting
fw <- function(y, z){
  # Change the flf formula so that all parameters are held in the vector, x
  ## Such as: x = c(x0, vf, k, n) (also above)
  ## This just needs to happen for Nelder-Mead to work, idk...
  neldermead_flf <- function(x){
    mean(abs(y - (x[2])/(1+exp(-x[3]*(z-x[1])))^(1/x[4])))
  }
  neldermead_flf
}
a <- fw(posVel$V2, posVel$V1) # need to change this with each new rawData
x0 = c(ix0, ivf, ik, iN) # only have ix0 values if run the processMasterFolder() from another script
b <- fminsearch(a, x0)
# b$optbase$xopt = fitted parameters location in the output
velFit_nm = flf(posVel$V1, b$optbase$xopt[[1,1]], b$optbase$xopt[[2,1]], b$optbase$xopt[[3,1]], b$optbase$xopt[[4,1]])
dataFit_nm = data.frame(posVel$V1, velFit_nm)

####################################################################################
# nls() Method
####################################################################################
m <- fflf(posVel$V2, posVel$V1, ix0, ivf, ik, iN)
velFit_nls <- flf(posVel$V1, m$x0, m$vf, m$k, m$n)
dataFit_nls = data.frame(posVel$V1, velFit_nls) %>%
  rename(pos = posVel.V1,
         fittedVel = velFit_nls)

####################################################################################
# Now Plot!
####################################################################################
ggplot(data = posVel, aes(x = posVel$V1)) +
  geom_point(aes(y = posVel$V2)) +
  geom_point(data = dataFit_nls, aes(y = fittedVel), color = "blue") + #nls() method
  geom_point(data = dataFit_nm, aes(y = velFit_nm), color = "red") + #neldermead/fminsearch() method
  ggtitle("Comparing Fitting Methods on Velocity Curve (red = Nelder-Mead, blue = nls)") +
  xlab("Position along Root (px)") +
  ylim(0, 5) +
  ylab("Velocity along Root (px/frame") +
  theme_bw()

####################################################################################

##########################################################################################################
# Below is code taken from flf.R, and in a conglomeration of trying to fit nls() and neldermead(). 
##The code above is a more concise version of this!
### This was done on 11.18.2020
##########################################################################################################

# Here on 11.17.2020, nls() fitting works well,
## Now need to add this to flf.R system to fit all RILs
# Hard coded x and y values:
# y = c(5, 42, 6, 2, 0)
# x = c(1, 2, 3, 4, 5)
# f <- y ~ m*x + b
# nls(f, start=list(m = 1, b = 1))

# Soft-coded x and y values, with real data:
# posVel = read.csv("~/Desktop/RILPop_tmp/RIL1/RIL1--RIL1_001_1.5--/rawData.csv", header = FALSE)
#
# # Best starting values for RIL1_005, using RIL1_001 parameters
# m <- fflf(posVel$V2, posVel$V1, ix0 = 350, ivf = 1.2, ik = 0.008, iN = 0.6)
# s <- summary(m)
# s$[["x0", "Estimate"]] # provides x0's value for that rawData.csv
# # Best starting values for RIL1_001 with guessing parameters
# m <- fflf(posVel$V2, posVel$V1, ix0 = 600, ivf = 1.5, ik = 0.001, iN = 1)

# data = data.frame(x, y)
# dataFit = data.frame(x, predict(m)) %>%
#       rename(pos = x,
#       fittedVel = predict.m.)
#
# ggplot(data = data, aes(x = x)) +
#   geom_point(aes(y = y)) +
#   geom_point(data = dataFit, aes(y = fittedVel), color = "blue") +
#   geom_point(data = dataFitY, aes(y = fitY), color = "red") +
#   ggtitle("Comparing Fitting Methods on Velocity Curve (red = Nelder-Mead, blue = nls)") +
#   xlab("Position along Root (px)") +
#   ylab("Velocity along Root (px/frame") +
#   theme_bw()


# Here on 11.18.2020
# Trying out neldermead function & loaded in neldermead package
# neldermead(x0, fn, lower = NULL, upper = NULL, nl.info = FALSE,
#            control = list(), ...)
# # where x0 = starting point, fn = function to minimize
# x <- c(x,x0,vf,k,n)
# fw <- function(y, z){
#   # Change the flf formula so that all parameters are held in the vector, x
#   ## Such as: x = c(x0, vf, k, n) (also above)
#   ## This just needs to happen for Nelder-Mead to work, idk...
#   neldermead_flf <- function(x){
#     mean(abs(y - (x[2])/(1+exp(-x[3]*(z-x[1])))^(1/x[4])))
#   }
#   neldermead_flf
# }
#
# x0 = c(ix0, ivf, ik, iN)
# neldermead(x0, neldermead_flf)



