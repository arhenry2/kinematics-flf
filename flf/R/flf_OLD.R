##### R Package for flf & plotting graphs w/ fit curve for imported data #####

####################################################################################
#Flexible Logistic Function: Function from Morris & Silk 1992 equation (Equation 8)
####################################################################################
flf <- function(x,x0,vf,k,n){
  (vf)/(1+exp(-k*(x-x0)))^(1/n)
}
####################################################################################

####################################################################################
#Flexible Logistic Function: Function from Morris & Silk 1992 equation (Equation 8)
####################################################################################
REGR <- function(x, x0, vf, k, n, scale){
  v = flf(x, x0, vf, k, n)
  (scale * k * v * ((vf^n) - (v^n)) ) / (n * (vf^n))
}
####################################################################################

####################################################################################
# SUMMARY: This function will find the folders in the masterPath
#           variable (string) and will "process" each folder
##########################################
# INPUTS: masterPath : string : the full path to the folder
#                     that will contain genotype folders
#                     such that each genotype folder will contain
#                     trial folders whereby each trial folder will
#                     contain a file by the name of rawData.csv
#                     After each folder is processed, this function will
#                     pairwise compare the results from the fits
##########################################
# OUTPUTS: obtained from flf_pairWise
##########################################
processMasterFolder <- function(masterPath){
  # this will find the folder names in the masterPath
  cdir <- dir(masterPath, NULL, FALSE, TRUE, FALSE, FALSE, TRUE);
  tmpR = vector('list')
  # for each folder - perform fit
  for (e in 1:length(cdir)) {
    grep(cdir[e], '/')
    # call processFolder
    resultT = processFolder(cdir[e])
    # store the result in a list of results
    tmpR[[e]] <- resultT
  }
  # perform pairwise compare between each fit result
  #processMasterFolder <- flf_pairWise(cdir,tmpR)
  tmpR
}
####################################################################################

####################################################################################
# SUMMARY: Loops through rawData.csv files from "pathName" (path to
#           folder w/ genotype replicates folders, which each
#           contain one rawData.csv file). Then creates an empty
#           dataframe to label parameters (coeffs). Creates empty
#           vector for the rawData values. for loop fits rawData.csv
#           via flf_fit_fromFile, labels K as coeffs, appends coeffs
#           from rawData.csv to parameter dataframe (as a row) to a
#           previously fitted rawData.csv, returns rawData in a list
#           to each parameter row.
##########################################
# INPUTS: pathName : full pathname to a folder that contains
#                     folders of analyzed (PhytoMorph) genotype
#                     replicates, and each of these folders have
#                     one rawData.csv file.
##########################################
# OUTPUTS: R: a summary table that contains a vector of lists:
#              first list - parameters of that rawData.csv file,
#              second list - pos values from rawData.csv,
#              third list - vel values from rawData.csv.
#            This collection of lists is produced from each
#            rawData.csv file run
##########################################
processFolder <- function(pathName){
  cdir <- dir(pathName, 'rawData.csv', FALSE, TRUE, TRUE);
  ktmp <- data.frame(fileName = c(), x0 = c(), vf = c(), k = c(), n = c())
  rawData <- c()
  for (e in 1:length(cdir)){
    print(cdir[e])
    D <- flf_fit_fromFile(cdir[e])
    K <- coef(D$coeffs)
    tmp <- data.frame(fileName = cdir[e], x0 = K[1], vf = K[2], k = K[3], n = K[4], mle = D$mle)
    ktmp <- rbind(ktmp,tmp)
    tmpRawList <- list(fileName = cdir[e], pos = D$pos, vel = D$vel)
    rawData[[e]] <- tmpRawList
  }
  R <- list(summary_table = ktmp, rawData = rawData)
}


####################################################################################
# SUMMARY: This function will load one rawData.csv file's pos and vel data
#           into a list. Assigns the first column as "pos" in R, and the 2nd
#           column as "vel"
##########################################
# INPUTS: fileName : csv file : the full path to the folder
#                     that contains one rawData.csv file
##########################################
# OUTPUTS: List containing vectors $pos and $vel from one rawData.csv file
##########################################
flfLoaderfromFile <- function(fileName){
  rawData <- read.csv(fileName, header = FALSE)
  pos <- rawData$V1
  vel <- rawData$V2
  ret <- list("pos" = pos, "vel" = vel)
  ret
}
####################################################################################


####################################################################################
# SUMMARY: This function will fit pos & vel values from csv file into flf, then
#           calculate the maximum likelihood estimation of initial parameters
#           to produce fitted parameters (parameters later labeled in flf_fit_fromFile).
##########################################
# INPUTS: inputList : csv file(s) : the full pathname to folder that contains
#                     rawData.csv file. Log likelihood of residuals measures
#                     parameters that are estimated via maximum likelihood estimation.
##########################################
# OUTPUTS: Fitted parameters (x0i, vfi, ki, and ni) from mle function
##########################################
fit_flf <- function(inputList){
  library(stats4)
  per = 0.0 # number we're using to define half width of the uniform distributions used for restarts (used in vfH, kH, nH) [Note: used as st dev for rnorm() case]
  pos <- inputList$pos
  vel <- inputList$vel
  svel <- sort(vel, decreasing = TRUE) #sort vf from highest to lowest values
  vfi <- mean(svel[1:100]) #grab highest 100 vf values
  vfH <- vfi*per # find half of vfi
  v0i <- mean(vel) #
  x0i <- mean(pos) #
  x0H <- x0i*per
  ni <- 1.2 #
  nH <- ni*per # find half of ni
  ki <- .008 #
  kH <- ki*per # find half of ki
  maxitr = 100 # added 11.19.2019 to see if number of iterations through mle changes the estimated parameter values (default = 100, had at 500 on 11.19.2019)
  reltol = 1e-20 # added 11.19.2019 to see if difference between values at interations changes the estimated parameter values (default = 1e-8, had at 1e-20 on 11.19.2019)
  reN = 5
  
  ll <- function(x0,vf,k,n){
    velp <- flf(pos,x0,vf,k,n)
    res <- vel-velp
    -sum(log(dnorm(res,0,1))) # old way before 3.2.2020
    # -sum(log(dnorm(res,0,0.1))) # Tried new one on 3.2.2020
    # -mean(log(dnorm(res,0,1)))
  }
  curMax = 0
  curMaxret = 0
  for (i in 1:reN){
    x0ii = rnorm(1, mean = x0i, sd = x0H) # define new range for x0i (add noise for where computer will start looking for max likelihood estimate for parameter x0)
    vfii = rnorm(1, mean = vfi, sd = vfH) # define new range for vfi (add noise for where computer will start looking for max likelihood estimate for parameter vf)
    kii = rnorm(1, mean = ki, sd = kH) # define new range for ki (add noise for where computer will start looking for max likelihood estimate for parameter k)
    nii = rnorm(1, mean = ni, sd = nH) # define new range for ni (add noise for where computer will start looking for max likelihood estimate for parameter n)
    start = list(x0 = x0ii, vf = vfi, k = ki, n = ni)
    ret <- mle(ll, start, method = "BFGS", control = list(maxit = maxitr, reltol = reltol))
    curV <- ll(coef(ret)[[1]], coef(ret)[[2]], coef(ret)[[3]], coef(ret)[[4]])
    if (curV > curMax){
      curMax = curV
      curMaxret = ret
    }
    
  }
  
  # ret <- mle(ll, start = list(x0 = x0i,vf=vfi, k=ki, n=ni), method = "BFGS", control = list(maxit = maxitr, reltol = reltol))
  # ll(coef(ret)[[1]], coef(ret)[[2]], coef(ret)[[3]], coef(ret)[[4]])
  
  # mle(ll, start = list(x0 = x0i,vf=vfi, k=ki, n=ni)) # Original mle function w/out testing parameters for outliers
  # mle(ll, start = list(x0 = x0i,vf=vfi, k=ki, n=ni), method = "BFGS", control = list(maxit = maxitr, reltol = reltol)) # New mle function that tests for outliers
  curMaxret = mle(ll, start = list(x0 = x0i,vf=vfi, k=ki, n=ni), method = "BFGS", control = list(maxit = maxitr, reltol = reltol)) # New mle function that tests for outliers
  
  list(coeffs = curMaxret, mle = curMax)
}
####################################################################################


####################################################################################
# SUMMARY: This function will load calculated fitted parameters from fit_flf.
##########################################
# INPUTS: fileName : csv file : the full path to the folder
#                     that contains rawData.csv file(s)
#                     The pos and vel values will be placed
#                     in separate vectors within the list
##########################################
# OUTPUTS: List containing vectors $pos and $vel from one rawData.csv file
##########################################
#takes file imported from analysis script, fits parameters to data
#"list" = returns data in sections: "$pos", "$vel", & "$coeffs"
flf_fit_fromFile <- function(fileName){
  outList <- flfLoaderfromFile(fileName)
  coeffs <- fit_flf(outList)
  
  list("pos" = outList$pos,"vel" = outList$vel,"coeffs" = coeffs$coeffs, "mle" = coeffs$mle)
}

####################################################################################
# SUMMARY: Takes data from flfDataStructure, calculates mean parameters and
#           standard deviation(!) barriers for each replicates. Each mean
#           parameter and st dev boundaries are plotted on the same graph,
#           one graph per replicate.
##########################################
# INPUTS: flfDataStructure : summary_table : the data structure is a list
#           of lists, each sublist containing the parameters for each
#           replicate and the rawdata for that replicate.
##########################################
# OUTPUTS: A dataframe that gives the average of each parameter for each
#           genotype.
##########################################


plot_flfDataStructure <- function(flfDataStructure, pos){
  
  # Store values in matrix of 4 rows, N columns [ N = # of replicates]
  sz     <- length(flfDataStructure$summary_table$x0)
  ku     <- matrix(list(), nrow=4, ncol=sz)
  ku[1,] <- flfDataStructure$summary_table$x0
  ku[2,] <- flfDataStructure$summary_table$vf
  ku[3,] <- flfDataStructure$summary_table$k
  ku[4,] <- flfDataStructure$summary_table$n
  
  ku[[1]] <- mean(flfDataStructure$summary_table$x0)
  ku[[2]] <- mean(flfDataStructure$summary_table$vf)
  ku[[3]] <- mean(flfDataStructure$summary_table$k)
  ku[[4]] <- mean(flfDataStructure$summary_table$n)
  
  # ku[[5]] <- sd(flfDataStructure$summary_table$x0)
  # ku[[6]] <- sd(flfDataStructure$summary_table$vf)
  # ku[[7]] <- sd(flfDataStructure$summary_table$k)
  # ku[[8]] <- sd(flfDataStructure$summary_table$n)
  
  return(ku)
}


###############################################################
################ Fitting curves to parameters #################
###############################################################

# This was written 2018 (when I was a wee babe) and for analyzing
## Cvi, Ler, Col, and pin2, abcb4, and abcb4pin2 mutants
# To change dataStructure from res[]:
# genotype <- #
# flfDataStructure <- res[[genotype]]
# where: 1 = abcb4, 2 = abcb4pin2, 3 = Col, 4 = Cvi, 5 = Ler, 6 = pin2, 7 = RIL1, 8 = RIL2, 9 = RIL3, and 10 = RIL4
## Number is the location in res[[]] (created in use_flf.R)

# genotype <- 1
# flfDataStructure <- res[[genotype]]
evaluateVelFits <- function(flfDataStructure, pos){
  
  # Store values in matrix of 4 rows, N columns [ N = # of replicates]
  sz     <- length(flfDataStructure$summary_table$x0)
  ku     <- matrix(list(), nrow=4, ncol=sz)
  ku[1,] <- flfDataStructure$summary_table$x0
  ku[2,] <- flfDataStructure$summary_table$vf
  ku[3,] <- flfDataStructure$summary_table$k
  ku[4,] <- flfDataStructure$summary_table$n
  
  # Store pos at 1st row then loop through replicates and get fitted values
  Gout <- matrix(list(), nrow = length(pos), ncol = sz+1)
  Gout[,1] <- pos
  for (n in 1:sz) {
    Gout[,n+1] <- flf(pos, ku[[1,n]], ku[[2,n]], ku[[3,n]], ku[[4,n]])
  }
  return(Gout)
}


evaluateREGRFits <- function(flfDataStructure, pos, scale){
  
  # Store values in matrix of 4 rows, N columns [ N = # of replicates]
  sz     <- length(flfDataStructure$summary_table$x0)
  ku     <- matrix(list(), nrow = 4, ncol = sz)
  ku[1,] <- flfDataStructure$summary_table$x0
  ku[2,] <- flfDataStructure$summary_table$vf
  ku[3,] <- flfDataStructure$summary_table$k
  ku[4,] <- flfDataStructure$summary_table$n
  
  # Store pos at 1st row then loop through replicates and get fitted values
  Gout <- matrix(list(), nrow = length(pos), ncol = sz+1)
  Gout[,1] <- pos
  for (n in 1:sz) {
    Gout[,n+1] <- REGR(pos, ku[[1,n]], ku[[2,n]], ku[[3,n]], ku[[4,n]], scale)
  }
  return(Gout)
}

# plots pos v vel w/ curve from parameters given from flf_fit_fromFile
#ypre: uses coefficients from mle as the fixed values for the curve on pos v vel
flf_plot <- function(pos, vel, coeffs){
  k <- coef(coeffs)
  spos <- sort(pos)
  ypre <- flf(spos, k[1], k[2], k[3], k[4])
  plot(pos, vel, col = "blue")
  flf_plotcurve(pos, coeffs)
}

# note this will be broken because k is passed in now
flf_plotcurve <- function(pos,k,HOLD,RANGY){
  #k <- coef(coeffs)
  spos <- sort(pos)
  ypre <- flf(spos,k[1],k[2],k[3],k[4])
  plot(pos,ypre,col="blue",new=HOLD,yim=RANGY)
  lines(spos,ypre)
}

################################################################################################
########################### Simulation Functions Written ~12.10.2018  ##########################
################################################################################################
### Simulates data using...?
sim <- function(x,sx,sv,x0,vf,k,n){
  vel <-flf(x,x0,vf,k,n)
  nx <- rnorm(length(x),0,sx) #rnorm gives numbers (mean of 0) to add/subtract (displacement) from real x values to make noise
  nv <- rnorm(length(x),0,sv) # " w/ y values
  vel <- vel+nv
  pos <- x+nx
  plot(vel~pos)
  list("V0"=x,"V1"=vel)
}

### Creates simulations for the parameters: x0, vf, k, and n ###
# p[x] referneces parameters in the array, 1 = x0, 2 = v0, 3 = k, 4 = n.
paraSim <- function(minPara, maxPara){
  n = 1
  p = vector('numeric')
  p[1] <- runif(n, min = minPara[1], max = maxPara[1])
  p[2] <- runif(n, min = minPara[2], max = maxPara[2])
  p[3] <- runif(n, min = minPara[3], max = maxPara[3])
  p[4] <- runif(n, min = minPara[4], max = maxPara[4])
  p
}

multiSim <- function(N,oPath,x,sx,sv,minPara,maxPara){
  for(i in 1:N){
    parameterValues <- paraSim(minPara,maxPara)
    tmpData <- sim(x,sx,sv,parameterValues[1],parameterValues[2],parameterValues[3],parameterValues[4])
    # Make main & sub directory (warning given if exists)
    fileName <- paste(oPath,as.character(i),sep="")
    fileName <- paste(fileName, "/", sep="")
    dir.create(fileName, recursive = 1)
    fileName <- paste(fileName, "rawData.csv", sep="")
    write.table(tmpData, file = fileName, row.names = FALSE, col.names=FALSE, sep=",")
  }
}

#copied flf_plotcurve from above, trying to plot avg parameters for Cvi
flf_plotavgcurve <- function(pos,coeffs){
  k <- coef(coeffs)
  spos <- sort(pos)
  ypre <- flf(spos,x0,vf,k,n)
  plot(pos,ypre,col="green")
  lines(spos,ypre)
}


################################################################################################
############################# T-test Parameters Written Awhile Back ############################
################################################################################################

#Function to t.test k and l parameters (Cvi v Ler for example)
compareResult <- function(k,l){
  parameterComparex0 <- t.test(k$x0,l$x0)
  parameterComparevf <- t.test(k$vf,l$vf)
  parameterComparek <- t.test(k$k,l$k)
  parameterComparen <- t.test(k$n,l$n)
  compareResultresult <- list(x0=parameterComparex0,vf=parameterComparevf,k=parameterComparek,n=parameterComparen)
  
}

#compares parameter output from file 1 v file 2, file 1 v file 3, etc through all files (files not named yet in output)
flf_pairWise <- function(cdir,tmpR){
  results = list();
  cnt = 1;
  for (i in 1:(length(tmpR)-1)){
    for (j in (i+1):length(tmpR)) {
      q1 <- basename(cdir[i])
      q2 <- basename(cdir[j])
      Q <- paste(q1,q2,sep='vs')
      results[[cnt]] <- list(compare=Q,results=compareResult(tmpR[[i]],tmpR[[j]]))
      cnt <- cnt + 1
    }
  }
  results
}
