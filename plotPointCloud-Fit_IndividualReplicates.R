

# Load in full res structure
load("~/rildata/rildata_realfullDataset/2020_11_19_res_RILPop_full.RData")

# Load in summary table for all RILs & all replicates
AllRILs_SummaryParameters = read.csv("~/rildata/rildata_realfullDataset/2020_11_20_AllRILs_SummaryParameters.csv")

# Density plot for mle values
ggplot(AllRILs_SummaryParameters, aes(x = mle)) +
  geom_density() +
  ggtitle("mle Density Plot for all RILs and Replicates")

# Select replicates with highest mean residual
highMleScores = AllRILs_SummaryParameters %>%
  filter(mle > 0.95)
lowMleScores = AllRILs_SummaryParameters %>%
  filter(mle < 0.06)

mleScores = AllRILs_SummaryParameters %>%
  filter(mle > 0.38)

# Load in chosen RIL & replicate, need to manually change w/ each new replicate
posVel = read.csv("/Users/ashleyhenry/Desktop/RILPop_tmp/RIL83_3/RILpop--RIL83_3--RIL83_001_4-- copy/rawData.csv", header = FALSE)
posVel = posVel %>%
  rename(
    pos = V1,
    vel = V2
  )

# Create fitted vel values with parameters in summary table
## Note: [1] = RIL42_003_6
##       [2] = RIL83_3--RIL83_001_4
##       [3] = RIL86_2--RIL86_003_4
# AllRILs_SummaryParameters$fileName[645] # Check replicate you're fitting!
# fittedVel = flf(posVel$pos, AllRILs_SummaryParameters$x0[645],	
#                 AllRILs_SummaryParameters$vf[645],	
#                 AllRILs_SummaryParameters$k[645],	
#                 AllRILs_SummaryParameters$n[645]) # fitted vel values for ground truth parameters

# Same as above, but with res:
res[[2]]$summary_table[2]$fileName
fittedVel = flf(posVel$pos, 
                res[[3]]$summary_table[2]$x0,	
                res[[3]]$summary_table[3]$vf,	
                res[[3]]$summary_table[4]$k,	
                res[[3]]$summary_table[5]$n)
fittedposVel = cbind(posVel$pos, fittedVel) # makes matrix of pos & fitted vel values for g truth parameters
fittedposVel = data.frame(fittedposVel) # forces the matrix to be a dataframe


# Lastly, plot point cloud & fitted line to check fit
ggplot(data = posVel, aes(x = posVel$pos)) + 
  geom_point(aes(y = vel)) +
  geom_point(data = fittedposVel, aes(y = fittedVel), color = "blue") +
  ggtitle("Velocity Fit for", res[[3]]$summary_table[1]$fileName) +
  xlab("Position from Root Tip (px)") +
  ylab("Velocity from Root Tip (px/frame)") +
  xlim(0,1500) +
  ylim(-1,5) +
  theme_bw()

