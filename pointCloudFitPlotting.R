# Plotting x0 values that I got v Nathan got to see if there's a correlation
## Update: there is not a correlation
x0Comparison = read.csv("~/Desktop/AshleyvNathan_x0.csv")
plot(x0Comparison$NATHAN_x0, x0Comparison$ASHLEY_x0)
mean(x0Comparison$NATHAN_x0, na.rm = TRUE)
mean(x0Comparison$Ashley_x0)


############################################################
####Script to plot flf fit over point cloud 10-30-2020 #####
############################################################

# Load in res structure with the full RIL set
load("~/rildata/rildata_fullDataset/res_RILs_full.RData")

# This res data doesn't include pos and vel values... Too big of a dataset?? Need to fix this!!!
## Just use rawData.csv from the folder?
## But the flf.R should already do this, no reason to do it twice...!

# only saves parameters, need to read.csv the rawData.csv from that folder
posVel = read.csv("~/Desktop/RILPop_tmp/RIL140_3/RIL140_3--RIL140_003_3.5--/rawData.csv")
posVel = posVel %>%
  rename(
    pos = X10.082,
    vel = X.0
  )

parameters = res[[1]]$summary_table

# Copied from AnalyzeKinematicAnalysisToolOutput.R, edited on 2020_10_27
## Creates flf ~fit~ line data, use the vel from these to plot the fit line
for (i in 1:length(res)) {
  resultVel <- evaluateVelFits(res[[i]], posVel$pos)
  matrixEnd = ncol(resultVel) - 1
  for (j in 2:matrixEnd) {
    resultVel[,j] <- as.numeric(resultVel[,j]) } # conversion from px/frame to mm/hr for all velocity value columns
  resultVel[,1] = as.numeric(resultVel[,1]) # conversion from px to mm for the one position value column
  write.csv(resultVel, paste('/Users/ashleyhenry/Desktop/PosEvaluatedVelRIL140_003_3.5', as.numeric(i), '.csv', sep = ""), row.names = F)
  plot(resultVel[,1], resultVel[,2])
  
  ######################################
  
  resultREGR <- evaluateREGRFits(res[[i]], posVel$pos, 120*100)
  matrixEnd = ncol(resultREGR) - 1
  
  # for (k in 2:matrixEnd) {
  # resultREGR[,j] <- as.numeric(resultREGR[,j]) * 120} # DOESN'T WORK: conversion from px/frame*frame to px/frame*hr (or %/hr) for all velocity value columns
  
  resultREGR[,1] = as.numeric(resultREGR[,1]) / 1463 # WORKS: conversion from px to mm for the one position value column
  write.csv(resultREGR, paste('/Users/ashleyhenry/Desktop/PosEvaluatedREGR', as.numeric(i), '.csv', sep = ""), row.names = F)
  plot(resultREGR[,1], resultREGR[,2])
}


data = read.csv("~/Desktop/PosEvaluatedVelRIL140_003_3.51.csv")
data = data %>%
  rename(
    pos = V1,
    fittedVel = V2
  )

# Plot point cloud w/ overlaying fitted velocity curve line
ggplot(posVel, aes(x=posVel$pos)) + 
  geom_point(aes(y = posVel$vel, color = "red")) +
  geom_line(aes(y = data$fittedVel))






