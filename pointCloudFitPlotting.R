# Plotting x0 values that I got v Nathan got to see if there's a correlation
## Update: there is not a correlation
x0Comparison = read.csv("~/Desktop/AshleyvNathan_x0.csv")
plot(x0Comparison$NATHAN_x0, x0Comparison$ASHLEY_x0)
mean(x0Comparison$NATHAN_x0, na.rm = TRUE)
mean(x0Comparison$Ashley_x0)


############################################################
####Script to plot flf fit over point cloud 10-30-2020 #####
############################################################
### Brain dump for what I want this script to do: ###
# co-plot rawData point clouds and fitted flf line (velocity plot, not REGR)
# use rawData part of res for point cloud
## or take rawData from res, with filename as the column names with pos & vel
# use code from AnalyzeKinematicAnalysisToolOutput.R to make fitted line out of given parameters in res
## have this save in the environment to use, not my Desktop
## has this save into res structure?
# be able to write in how many RILs or files I want plotted, so this needs to be a ~function~!
# have this highlight bad fits (R-squared, other outlier tests using data and fitted line)
## look up how to find a bad fit...
## Logistic functions good fit resource: https://statisticalhorizons.com/wp-content/uploads/GOFForLogisticRegression-Paper.pdf
## General model good-fit resource: https://www.itl.nist.gov/div898/handbook/pmd/section4/pmd44.htm < likely the better resource for flf

functionName = function(argument, argument, argument){
  # Grab rawData for replicate(s) - rawData.csv in replicate's folder
  # Grab parameters for replicate(s) - in summary table
  # Calculate fitted line from parameters for replicate(s) - calc from summary table^
  # Plot rawData pos & vel
  #   plot fitted line

  

}
#####################################################################################
################################# Here on 11.3.2020! ################################
#####################################################################################
# Goal: take out rawData from res structure & save into own dataframe to be used in plotting function above
tmp2 <- as.data.frame(res[[1]]$rawData[[i]]) # This makes a dataframe with pos and vel columns! But filename is it's own column...
# Now use this ^ in a for loop

DataframeToFill = matrix(1, nrow(length(res[[1]]$rawData))) # Error here with nrow()
outnames = matrix(1, nrow(genotype)) # example code without errors

emptyDataframe = data.frame()
for (i in 1:length(res[[1]]$rawData)){
  tmp2 <- as.data.frame(res[[1]]$rawData[[i]]) # This makes a dataframe with pos and vel columns! But filename is it's own column...
  
}


# - Plot flf fit on these point clouds, because they're quirky
## - finished this on 11.2.2020 & saved plots to Desktop
  # RIL1--RIL1_009_2.5--
  # RIL1--RIL1_010_2--
  # RILpop--RIL1--RIL1_003_2.5--

# Load in res structure with the full RIL set
# load("~/rildata/rildata_fullDataset/2020_10_30_res_RILs_full.RData")
load("~/rildata/rildata_fullDataset/2020_11_2_res_RIL1.RData")

# only saves parameters, need to read.csv the rawData.csv from that folder
posVel = read.csv("~/Desktop/RILPop_full/RIL1_1/RIL1--RIL1_009_2.5--/rawData.csv")
posVel = posVel %>%
  rename(
    pos = X10.046,
    vel = X0
  )


# Copied from AnalyzeKinematicAnalysisToolOutput.R, edited on 2020_10_27
## Creates flf ~fit~ line data, use the vel from these to plot the fit line

# Brain Dump:
# - have for loop take fileName from rawData & use that for column name for each replicate - this has been a bitch...

pos <- seq(0, 2000, 1)

for (i in 1:length(res)) {
  resultVel <- evaluateVelFits(res[[i]], pos) # uses parameters to calc fitted velocity curve, converts px/frame to mm/hr in that function
  matrixEnd = ncol(resultVel) - 1
  for (j in 2:matrixEnd) {
    resultVel[,j] <- as.numeric(resultVel[,j]) # fitted vel values
    }
    
  resultVel[,1] = as.numeric(resultVel[,1]) # pos values
  # resultVel = data.frame(resultVel)
  # colnames(resultVel) <- paste(res[[i]]$rawData[[j]]$fileName)
  write.csv(resultVel, paste('/Users/ashleyhenry/Desktop/PosFittedVel_RIL1_', as.numeric(i), '.csv', sep = ""), row.names = F)
  plot(resultVel[,1], resultVel[,2])
  
  # The code to do the same as above for a REGR curve is located in: AnalyzeKinematicAnalysisToolOutput.R
}


# Fitted Velocity Values for Fitted Line in Plot
data1 = read.csv("~/Desktop/PosFittedVel_RIL1_1.csv")
data2 = read.csv("~/Desktop/PosFittedVel_RIL1_2.csv")
data1 = data1 %>%
  mutate(pos = pos/1500) %>%
  mutate(vel = V10*0.08202)

data1 %>%
  for (i in 1:length(res)) {
    matrixEnd = ncol(data1) - 1
    for (j in 2:matrixEnd) {
      colnames(data1$V[i]) <- paste(res[[i]]$rawData[[j]]$fileName)
      rename(
        pos = V1)
    }
  }


# data2 = read.csv("~/Desktop/PosFittedVel_RIL1_2.csv")

RIL1_fittedVelData = cbind(data1, data2)


# Plot point cloud w/ overlaying fitted velocity curve line
ggplot(posVel, aes(x = pos)) + 
  geom_point(aes(y = vel)) +
  geom_point(data = data1, aes(x = V1, y = V10), color = "blue") +
  ggtitle("Velocity Fit for RIL1_009_2.5") +
  xlab("Position from Root Tip (px)") +
  ylab("Velocity from Root Tip (px/frame)") +
  xlim(0,1500) +
  theme_bw()



########################################################
# Plotting point cloud w/ flf fit line overlaying that #
########## Taken from RILSetAnalysisScript.R ###########
########################################################

#######################################################
# Saving flf fit for that RIL140_3 > RIL140_003_3.5
# parameters = res[[2]]$summary_table
parameters = RES[[1]][[106]]$summary_table

# parameters = res[[2]]$summary_table
# parameters = res[[3]]$summary_table
# parameters = res[[4]]$summary_table
# parameters = res[[5]]$summary_table
# parameters = res[[6]]$summary_table

# Averages the parameters, then fit velocty values to that pos and avg parameters
## Table that gives avg and sd for each paramters in that condition/genotype
parameter_summary = parameters %>%
  summarize(avg_x0 = mean(parameters$x0),
            avg_vf = mean(parameters$vf),
            avg_k = mean(parameters$k),
            avg_n = mean(parameters$n))

# Matrix that has position values in 1st column and fitted velocity values in 2nd column
avgVelCurve <- matrix(nrow = length(pos), ncol = 2)
avgVelCurve[,1] <- pos
avgVelCurve[,2] <- flf(pos, parameter_summary$avg_x0, parameter_summary$avg_vf, parameter_summary$avg_k, parameter_summary$avg_n)

write.csv(avgVelCurve, paste('/Users/ashleyhenry/Desktop/avgVelCurve_RIL140_3_003_3.5.csv'), row.names = F)



