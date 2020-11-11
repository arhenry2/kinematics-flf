# New R Script to fit my groung truth parameters
## Copied fit_flf() from flf.R, but took it out of a function to get my fitted parameters from the ground truth
## Having them in a folder using the masterPath & res loading isn't working, so I'll have to break down the steps to see what's wrong
fit_flf <- function(inputList){
  library(stats4)
  per = 0.1 # number we're using to define half width of the uniform distributions used for restarts (used in vfH, kH, nH) [Note: used as st dev for rnorm() case]
  pos <- gTruthPosVel$pos
  vel <- gTruthPosVel$vel
  svel <- sort(vel, decreasing = TRUE) #sort vf from highest to lowest values
  vfi <- mean(svel[1:100]) #grab highest 100 vf values
  print(c("mean of highest vf values", vfi))
  vfH <- vfi*per # find half of vfi
  v0i <- mean(vel) #
  print(c("mean of vel values", v0i))
  x0i <- mean(pos) #
  print(c("mean of pos values", x0i))
  x0H <- x0i*per
  ni <- 1.2 #
  nH <- ni*per # find half of ni
  ki <- .008 #
  kH <- ki*per # find half of ki
  # maxitr = 100 # added 11.19.2019 to see if number of iterations through mle changes the estimated parameter values (default = 100, had at 500 on 11.19.2019)
  maxitr = 10000
  # reltol = 1e-20 # added 11.19.2019 to see if difference between values at interations changes the estimated parameter values (default = 1e-8, had at 1e-20 on 11.19.2019)
  reltol = 1e-100
  reN = 5
  ll <- function(x0,vf,k,n){
    velp <- flf(pos,x0,vf,k,n)
    res <- vel-velp # calc's residuals
    # -sum(log(dnorm(res,0,1))) # old way before 3.2.2020
    -mean(log(dnorm(res, 0, 0.1))) # new way on 11.10.2020, gives probability of your residuals falling into a N dist (assuming N dist.)
    # -sum(log(dnorm(res,0,0.1))) # Tried new one on 3.2.2020
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
    start = list(x0 = x0ii, vf = vfi, k = ki, n = ni)
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
  curMaxret = mle(ll, start = list(x0 = x0i, vf=vfi, k=ki, n=ni), method = "BFGS", control = list(maxit = maxitr, reltol = reltol)) # New mle function that tests for outliers
  
  # list(coeffs = curMaxret, mle = curMax)
  coeffs <- list(coeffs = curMaxret, mle = curMax)
  }
  print(c("coeffs:", curMaxret))
  # return(curMaxret)
  coeffs
}

fitTryData <- fit_flf(gTruthPosVel)


# Makes gTruthPosVel with correct column names for flf_fit()
pos <- seq(0, 1000, 1) # pos values for ground truth parameters
# fitted vel values for ground truth parameters
## Ground truth parameters from RIL1_1/RIL1--RIL1_003_4--!!!
vel = flf(pos, 658.9043,	1.0581177,	0.014531372,	2.9281828) 
gTruthPosVel = cbind(pos, vel) # makes matrix of pos & fitted vel values for g truth parameters
gTruthPosVel = data.frame(gTruthPosVel) # forces the matrix to be a dataframe

velFit = flf(pos, fitTryData[[1]]@coef[[1]], fitTryData[[1]]@coef[[2]], fitTryData[[1]]@coef[[3]], fitTryData[[1]]@coef[[4]])
# plot(pos, velFit)

### Add fitted vel from fit_flf() above
gTruthPosVelFit = cbind(pos, vel, velFit)
gTruthPosVelFit = data.frame(gTruthPosVelFit)

ggplot(gTruthPosVelFit, aes(x = pos)) +
  geom_point(aes(y = vel), color = "blue") +
  geom_point(aes(y = velFit)) +
  theme_bw()



gTruth
# "gTruth" gives the following parameters:
#  x0           vf            k            n 
# 658.90430000   1.05811770   0.01453137   2.92818280
ret
# "ret" gives the following parameters:
#  x0           vf            k            n 
# 999.79974068   1.10537089   0.01387277   7.81947478 














  