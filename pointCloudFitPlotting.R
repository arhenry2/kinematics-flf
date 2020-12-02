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

rawDataSummaryTable = data.frame(matrix( # Error with nrow(), so wrote 3 rows: filename, pos, vel
                             ncol = length(res[[1]]$rawData[[1]]),
                             nrow = 0,
                             dimnames = list(NULL, c("filename", "pos", "vel"))))

for (i in 1:length(res[[1]]$rawData)) {
  # gtype = res[[i]]$summary_table
  # summaryTable = gList[[1]]
  # rawDataSummaryTable <- as.data.frame(res[[1]]$rawData[[i]])
  rawData <- as.data.frame(res[[1]]$rawData[[i]]) # filename, pos, & vel have columns, how to make filename part of the column name?!
  rawDataSummaryTable <- cbind(rawData, rawDataSummaryTable) # trying to bind rawData to summary table; won't work b/c dif # of rows
  # rawDataSummaryTable %>%
  #   # mutate(res[[1]]$rawData[[i]]$fileName_pos <- res[[1]]$rawData[[i]]$pos) # want the filename to be in the column name
}

for (i in 1:length(res[[1]]$rawData)){
  tmp2 <- as.data.frame(res[[1]]$rawData[[i]]) # This makes a dataframe with pos and vel columns! But filename is it's own column...
}

# - Plot flf fit on these point clouds, because they have low MLE values in the SummaryTable
## RILpop--RIL1--RIL1_004_3-- w/ MLE = 8897.891 - DONE
## RILpop--RIL135--RIL135_002_5--	w/ MLE = 8440.048 - DONE
## RIL79_2--RIL79_003_4.5--	w/ MLE = 9982.179		- DONE

# Load in res structure with the full RIL set
load("~/rildata/rildata_fullDataset/2020_11_12_res_RIL1.RData")

pos <- seq(0, 2000, 1)
flfDataStructure = res[[1]]
sz     <- length(flfDataStructure$summary_table$x0)
ku     <- matrix(list(), nrow = 4, ncol = sz)
ku[1,] <- flfDataStructure$summary_table$x0
ku[2,] <- flfDataStructure$summary_table$vf
ku[3,] <- flfDataStructure$summary_table$k
ku[4,] <- flfDataStructure$summary_table$n

# Store pos at 1st row then loop through replicates and get fitted values
Gout <- matrix(list(), nrow = length(pos), ncol = sz)
# Gout[,1] <- pos

# for (i in 1:length(res[[1]]$rawData)){
  # print(res[[1]]$rawData[[i]]$fileName) # works
  for (n in 1:sz) {
    Gout[,n] <- flf(pos, ku[[1,n]], ku[[2,n]], ku[[3,n]], ku[[4,n]])
    plot(res[[1]]$rawData[[n]]$pos, res[[1]]$rawData[[n]]$vel,
         main = paste(res[[1]]$rawData[[n]]$fileName),
         xlab = "Pos (px)",
         ylab = paste("Vel for ", n)) # works
    lines(pos, Gout[,n],
          col = "blue")
  }  
# }
###############################################################
######### Ground Truth Parameters | 10 November 2020 ##########
###############################################################
# Need to find make a ground truth & see if it stays that way!
## Pick parameters & make fitted line - done
## Add noise to this fitted line to make a point cloud
## See if you can get the same parameters from that point cloud
gTruth = c(x0 = 658.9043, vf = 1.0581177, k = 0.014531372, n = 2.9281828)
# ground truth parameters from RIL1_1/RIL1--RIL1_003_4--!!!
pos <- seq(0, 2000, 1) # pos values for ground truth parameters
fittedVel = flf(pos, 658.9043,	1.0581177,	0.014531372,	2.9281828) # fitted vel values for ground truth parameters
gTruthFitted = cbind(pos, fittedVel) # makes matrix of pos & fitted vel values for g truth parameters
gTruthFitted = data.frame(gTruthFitted) # forces the matrix to be a dataframe
plot(pos, fittedVel, main = "Ground Truth flf", xlab = "Position (px)", ylab = "Vel (px/frame)")
ggplot(gTruthFitted, aes(x = pos)) + # See if ggplot makes the fitted line look different, update: it doesn't
  geom_point(aes(y = fittedVel)) +
  ggtitle("GGplot Version") +
  theme_bw()

# Fit without noise:
write.csv(gTruthFitted, file = "~/Desktop/RILPop_tmp/RIL1_1/rawData.csv")
###############################################################
###############################################################
###############################################################


# only saves parameters, need to read.csv the rawData.csv from that folder
posVel = read.csv("/Users/ashleyhenry/Desktop/RILPop_tmp/RIL1_1/RIL1--RIL1_005_3--/rawData.csv", header = FALSE)
posVel = posVel %>%
  rename(
    pos = V1,
    vel = V2
  )


# Copied from AnalyzeKinematicAnalysisToolOutput.R, edited on 2020_10_27
## Creates flf ~fit~ line data, use the vel from these to plot the fit line

pos <- seq(0, 1500, 1)
for (i in 1:length(res)) {
  resultVel <- evaluateVelFits(res[[i]], pos) # uses parameters to calc fitted velocity curve
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

# Trying above code, but for one list of res:
## This works! Just need to change the res location for the replicates you're looking for! 2020-11-10
resultVel <- evaluateVelFits(res[[1]], pos)
matrixEnd = ncol(resultVel) - 1
for (j in 2:matrixEnd) {
  resultVel[,j] <- as.numeric(resultVel[,j]) # fitted vel values
}
resultVel[,1] = as.numeric(resultVel[,1]) # pos values
write.csv(resultVel, paste('/Users/ashleyhenry/Desktop/PosFittedVel_RIL1.csv', sep = ""), row.names = F)

# Fitted Velocity Values for Fitted Line in Plot
# data1 = read.csv("~/Desktop/PosFittedVel_RIL1_1.csv")
# data2 = read.csv("~/Desktop/PosFittedVel_RIL1_2.csv")
# data2 = data2 %>%
#   mutate(pos = pos/1500) %>%
#   mutate(vel = V5*0.08202)
# 
# data1 %>%
#   for (i in 1:length(res)) {
#     matrixEnd = ncol(data1) - 1
#     for (j in 2:matrixEnd) {
#       colnames(data1$V[i]) <- paste(res[[i]]$rawData[[j]]$fileName)
#       rename(
#         pos = V1)
#     }
#   }
# 
# RIL1_fittedVelData = cbind(data1, data2)

data = read.csv("/Users/ashleyhenry/Desktop/PosFittedVel_RIL1_1.csv")
# data = data %>%
#   mutate(pos = pos/1500) %>%
#   mutate(vel = V6*0.08202)

# Plot point cloud w/ overlaying fitted velocity curve line
ggplot(posVel, aes(x = pos)) + 
  geom_point(aes(y = vel)) +
  geom_point(data = data, aes(x = V1, y = V6), color = "blue") +
  ggtitle("Velocity Fit for RIL1_005_3") +
  xlab("Position from Root Tip (px)") +
  ylab("Velocity from Root Tip (px/frame)") +
  xlim(0,1500) +
  theme_bw()

## RILpop--RIL1--RIL1_004_3-- w/ MLE = 8897.891 - DONE
## RILpop--RIL135--RIL135_002_5--	w/ MLE = 8440.048 - DONE
## RIL79_2--RIL79_003_4.5--	w/ MLE = 9982.179	- DONE
## RILpop--RIL109--RIL109_003_3--	w/ MLE = 5913.303	
## RILpop--RIL96--RIL96_004_3.5--	w/ MLE = 9018.691
## RIL96_2--RIL96_005_4.5--	w/ MLE = 7298.717
## RILpop--RIL139--RIL139_001_3-- w/ MLE = 9377.972	
## RILpop--RIL144X--RIL144X_003_3--	w/ MLE = 9671.525	
## RILpop--RIL146_2--RIL146_004_4--	w/ MLE = 7064.311	
## RILpop--RIL41--RIL41_002_4--	w/ MLE = 6829.046


allRILsSummaryTable = read.csv("~/rildata/rildata_fullDataset/2020_10_21_AllRILs_SummaryParameters.csv")
# dim(allRILsSummaryTable) = 1607 x 8


allRILsSummaryTable_lowMLE = allRILsSummaryTable %>%
  filter(mle < 10000)
# dim(allRILsSummaryTable_lowMLE) = 156 x 8
# Random replicates I'm plotting to see if actually bad fits:

## RILpop--RIL1--RIL1_004_3-- w/ MLE = 8897.891 - DONE
## RILpop--RIL135--RIL135_002_5--	w/ MLE = 8440.048 - DONE
## RIL79_2--RIL79_003_4.5--	w/ MLE = 9982.179	- DONE
## RILpop--RIL109--RIL109_003_3--	w/ MLE = 5913.303	
## RILpop--RIL96--RIL96_004_3.5--	w/ MLE = 9018.691
## RIL96_2--RIL96_005_4.5--	w/ MLE = 7298.717
## RILpop--RIL139--RIL139_001_3-- w/ MLE = 9377.972	
## RILpop--RIL144X--RIL144X_003_3--	w/ MLE = 9671.525	
## RILpop--RIL146_2--RIL146_004_4--	w/ MLE = 7064.311	
## RILpop--RIL41--RIL41_002_4--	w/ MLE = 6829.046

# 156/1607 = 0.0970753
## The portion of fits that might be bad = ~10%







