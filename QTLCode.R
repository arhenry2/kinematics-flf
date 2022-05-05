# install.packages("data.table")
# install.packages("qtl2", repos="http://rqtl.org/qtl2cran")

library(data.table)
library(qtl2)
# library(qtl)

# Folder with all of the input files for this R script
setwd("~/Desktop/QTL_AnalysisMaterials/Ashley_QTL")

################# R/qtl on CvixLer RIL Population, Kinematic parameters VF, K, N ################
#/CvixLer    
#/out                         - output folder from analysis
# /CvixLer.working.json       - JSON which calls the objects dependents for rQTl2
# /CvixLer.RIL.MAP.csv        - MAP information, with SNP name, Chr, and Pos for each sitename
# /CvixLer.RIL.SNPs.csv       - SNPs in numeric coded format for R/qtl, 
# /CvixLer.RIL.AVGphe.csv     - Phenotypes
################################### Mapping #####################################################
# read json object which is calling the MAP, NUM, and Pheno within rQTLmeta
# CT <-read_cross2("CvixLer.working.json") # can use this if in correct directory
  CT <- read_cross2("~/Desktop/QTL_AnalysisMaterials/Ashley_QTL/CvixLer.working.json")
  CT   # summary; crosstype "riself" for 2-way RIL by selfing 
  # names(CT)
  # head(CT$geno)
  # head(CT$pheno)
# CT$geno is a list of 5 chromosomes, combine 5 chromosomes via cbind 
## do.call() allows cbind of multiple objects/vectors
  g1 <- do.call("cbind", CT$geno) 
# Pseudomarkers improve resolution, but need huge memory; step=0 does nothing
  map <- insert_pseudomarkers(CT$gmap, step = 2.5, stepwidth = "max")
  pr <- calc_genoprob(CT, map, error_prob = 0.002)
# For each individual at each position, find the genotype with the maximum marginal probability
## max marginal probability = the max prob of an event irrespective of the outcome of another variable (instead of conditional prob)
## I think it converts "pr" info into "1" or "NA"?
  m1 <- maxmarg(pr)
# Estimate the numbers of crossovers in each individual on each chromosome
  xo1 <- count_xo(m1) # produced dataframe with RILs vs Chr 1-5, filled in w/ zeros (no crossovers?)
# Calculate kinship matrix for individuals (accounting for relationship among indiv, i.e. including polygenic effect)
# loco = Leave One Chromosome Out, method (scan each chromosome using a kinship matrix that is calculated using data from all other chromosomes)
  k_loco <- calc_kinship(pr, "loco") # produced list of lists containing #s from 0.7-0.9 (what is this??)

## Perform genome scan by Haley-Knott regression
# Output of scan1() is matrix of LOD scores, positions x phenotypes
# Selecting phenotype columns & make sure to indicate the same columns in scan1perm() below
  out <- scan1(pr, CT$pheno[,c(1:4)], k_loco, cores = 4) ## Note: (1:4) for REGR Curve Descriptors, (1,2,8,9) for parameters

  summary(out)
  nperm = 10000
# Permutation testing
  operm <- scan1perm(pr[,1:5], CT$pheno[,c(1:4)], k_loco[1:5], addcovar = NULL,
                   Xcovar = NULL, intcovar = NULL, weights = NULL, reml = TRUE,
                   model = "normal", n_perm = nperm, perm_Xsp = FALSE, # changed n_perm=10000 from 1000 (5.9.2021)
                   perm_strata = NULL, chr_lengths = NULL, cores = 4)

  # summary(operm)
  # sig <- summary(operm, alpha = c(0.05, 0.01)) # Significant threshold, 95% & 99%  
  sig <- summary(operm, alpha = c(0.05)) # Significant threshold, just 95%

##### Finding LOD peaks

  peaks <- find_peaks(out, map, threshold = sig[1,1], drop = 1.5)
  # peaks_fullDataset <- find_peaks(out, map, threshold = sig[1,1], drop = 1.5)
  # peaks_noMatZone <- find_peaks(out, map, threshold = sig[1,1], drop = 1.5)
  # peaks_littleMatZone <- find_peaks(out, map, threshold = sig[1,1], drop = 1.5)
  
  peaks

  ### Testing Bayes peaks (if there are multiple peaks on one chromosome, like parameters k and n on chr 1)
  ## Saw this on Broman's website: https://kbroman.org/qtl2/assets/vignettes/user_guide.html#Finding_LOD_peaks,
  ### He says to use this with caution... BUT WHY?? WHAT IS BAYES?! WHAT IS PEAKDROP?!
  # Please Note: you have to change lodcolumn & chr manually to find the right peaks
  lod_int(out, map, lodcolumn = 2, chr = 4, peakdrop = 1.8, drop = 1.5) # Also works on lodcolumn = 2 (vf)

  ## Use this if want to normalize traits to have same threshold at 1! (added 4-19-2021)
  # for loop to divide each trait's data by that trait's sig value at 95%
  # Note: ONLY DO THIS ONCE! or it'll keep changing 'out'
  # for (i in 1:dim(out[,])[2]){
  #   out[,i] = out[,i]/sig[,i]
  # }
  
  
  # Trying out Interaction from Broman's FAQ page
  ## 7 March 2022
  # data(hyper)
  # hyper <- calc.genoprob(hyper, step=2.5, err=0.001)
  # pr <- calc.genoprob(CT, step = 2.5, err = 0.002)
  # 
  # 
  # qtl <- makeqtl(hyper, chr=15, pos=18, what="prob")
  # 
  # out.i <- addqtl(hyper, qtl=qtl, formula=y~Q1*Q2, method="hk")
  # out.a <- addqtl(hyper, qtl=qtl, formula=y~Q1+Q2, method="hk")
  # plot(out.i - out.a)
  # 
  # out2 <- scantwo(hyper, verbose=FALSE)
  # out3 <- scantwo(pr, CT$pheno[,c(5:8)], k_loco)
  # plot(out2, chr=c(1,4,6,7,15)) # pg. 217 of Broman's R/qtl manual
  # plot(out2, chr=c(1,4,6,7,15), lower="cond-int") # cleans up the plot to highlight pairs of QTL
  # plot(out2, chr=1, lower="cond-int", upper="cond-add")
  # 
  # # Trying it out on my data:
  # library(qtl)
  # CT_new <- read.cross("csvs", dir="~/Desktop/Mydata", genfile="ALFP_geno.csv",
  #                      phefile="2022-02-07_predGRfromWidth-maxREGR.csv", mapfile="ALFP_gmap.csv")
  # pr <- calc.genoprob(CT_new, map, step=2.5, err=0.001)
  # qtl <- makeqtl(CT, chr = 3, pos=2.613844, what="prob")
  # out.i <- addqtl(CT, qtl = qtl, formula=y~Q1*Q2, method="hk")
  # out.a <- addqtl(CT, qtl = qtl, formula=y~Q1+Q2, method="hk")
  # plot(out.i - out.a)
  
  
# Broman's for loop for plotting LOD scores
# Transform these with your data to easier plotting!)
  # color <- c("slateblue", "violetred", "green3")
  # par(mar=c(4.1, 4.1, 1.6, 1.1))
  # ymx <- max(maxlod(out), maxlod(out_pg), maxlod(out_pg_loco))
  # for(i in 1:2) {
  #   plot(out, map, lodcolumn=i, col=color[1], main=colnames(iron$pheno)[i],
  #        ylim=c(0, ymx*1.02))
  #   plot(out_pg, map, lodcolumn=i, col=color[2], add=TRUE)
  #   plot(out_pg_loco, map, lodcolumn=i, col=color[3], add=TRUE, lty=2)
  #   legend("topleft", lwd=2, col=color, c("H-K", "LMM", "LOCO"), bg="gray90", lty=c(1,1,2))
  # }
  
##### QTL Mapping Plots:
  ymx <- maxlod(out)
  plot(out, map, lodcolumn = 1, ylim = c(0, ymx*1.02), col = rgb(0, 0, 0, 0.5), main = "") # black
  plot(out, map, lodcolumn = 2, col = rgb(1, 0, 0, 0.5), add = TRUE) # red
  plot(out, map, lodcolumn = 3, col = rgb(1, 0.4, 0, 0.5), add = TRUE) # orange
  plot(out, map, lodcolumn = 4, col = rgb(0, 0.6, 1, 0.5), add = TRUE) # light blue
  # plot(out, map, lodcolumn = 5, col = rgb(0, 1, 0, 0.5), add = TRUE) # green
  # plot(out, map, lodcolumn = 6, col = rgb(0, 0, 0.8, 0.5), add = TRUE) # dark blue
  # plot(out, map, lodcolumn = 7, col = rgb(0.4, 0, 0.6, 0.5), add = TRUE) # violet
  # plot(out, map, lodcolumn = 8, col = rgb(0, 0, 0, 0.3), add = TRUE) # grey

  # MAKE ABLINE AT 1 AS THRESHOLD FOR ALL! (AFTER YOU NORMALIZE DATA BY SIG VALUES!) #added 4-19-2021
  # abline(h = 1, col = rgb(0, 0, 0, 0.5), lty = 1)
  abline(h = sig[,1], col = rgb(0, 0, 0, 0.5), lty = 2) # black
  abline(h = sig[,2], col = rgb(1, 0, 0, 0.5), lty = 2) # red
  abline(h = sig[,3], col = rgb(1, 0.4, 0, 0.5), lty = 2) # orange
  abline(h = sig[,4], col = rgb(0, 0.6, 1, 0.5), lty = 2) # light blue
  # abline(h = sig[,5], col = rgb(0, 1, 0, 0.5), lty = 2) # green
  # abline(h = sig[,6], col = rgb(0, 0, 0.8, 0.5), lty = 2) # dark blue
  # abline(h = sig[,7], col = rgb(0.4, 0, 0.6, 0.5), lty = 2) # violet
  # abline(h = sig[,8], col = rgb(0, 0, 0, 0.1), lty = 2) # grey
  legend("topright",
         lwd = 3,
         col = c(rgb(0, 0, 0, 0.5), # black
                 rgb(1, 0, 0, 0.5), # red
                 rgb(1, 0.4, 0, 0.5), # orange
                 rgb(0, 0.6, 1, 0.5)), # light blue
                 # rgb(0, 1, 0, 0.5), # green
                 # rgb(0, 0, 0.8, 0.5), # dark blue
                 # rgb(0.4, 0, 0.6, 0.5), # violet
                 # rgb(0, 0, 0, 0.3)), # grey
         colnames(out), bty = "n")

# plot(out, map[peaks[2,3]], lodcolumn = 2, xlim = c(peaks[2,6],peaks[2,7]),
                                      # ylim = c(0,7), main = paste0(colnames(out)[2]))
dev.off()

#####################################################################################################
# Create file with zeros and ones to show confidence intervals that are significant for each trait
## Added 2022-02-08 with Dr. Nathan Miller
# if LOD > sig[,#], then add 1
# if LOD < sig[,#], then add 0

# make a matrix to put ones and zeros:
ohYeah = matrix(out, nrow = nrow(out), ncol = ncol(out))

for (i in 1:length(sig)){
    for (j in 1:dim(out)[1]){
  
    if(out[j,i] >= sig[i]){
      # if that value is sig, add "one" in the matrix
      ohYeah[j,i] = 1
    } else {
      # if that value isn't sig, add "zero" in the matrix
      ohYeah[j,i] = 0
    }  
  }
}
# Change column names so you know what's what
colnames(ohYeah) = colnames(out)

write.csv(ohYeah, "~/Desktop/2022-02-14_sigPeaks_basis.csv")

#####################################################################################################
# Ashley creating csv file from "out" & "map"
out
map[1] # = chr1, etc. through chr 5
markerPositions = unlist(cbind(map[1], map[2], map[3], map[4], map[5]))
markerPositions = data.frame(markerPositions)
out = data.frame(out)
cm_Positions = cbind(markerPositions, out)
write.csv(cm_Positions, "~/Desktop/2022-02-14_cMPositions_FINAL.csv")

#####################################################################################################
# Load in geno data & see Landsberg erecta occurrences & pick all lines w/ Ler at Marker 1
## Note: This only extracts occurrences of L, it doesn't connect the phenotype at that marker location;
###     I ditched this idea 10-20-2020 when I found how to make pheno v geno plots (below: plot_pxg())
# library(tidyverse)
# geno = read.csv("~/Desktop/QTL_WuPracticewithSteveDeslauriers/Ashley_QTL/ALFP_geno.csv")
# geno1L = geno %>%
#   select(PVV4) %>%
#   filter(PVV4 == 'L')

#####################################################################################################
################################# Allele Effect at Specified Marker #################################
############################################ 2020-10-20 #############################################
# Finds if there are multiple peaks on one chr for each phenotype
# lod columns: 1:maxREGR, 2:vf, 3:positionMaxREGR, 4:growthZoneWidth
# This is a matrix, how can I add chr & phenotype to this & add to 'peaks' dataframe?
## Added this^ manually to peak data: 2021-05-09_peaks.xlsx
lod_int(out, map, lodcolumn = 4, chr = 1, peakdrop = 1.8, drop = 1.5) 
# for kinematic traits: lodcolumn = 2 & 4 for chr = 1

# maxmarg() grabs allele data by RIL at a specific marker location (given chromosome & position)
g <- maxmarg(pr, map, chr = 5, pos = 76.680571, return_char = TRUE, minprob = 0)
# Plots the phenotypes and genotype at a marker location
## bgcolor = 0 changes background color to white
plot_pxg(g, CT$pheno[,"overallGrowthRate"], ylab = "Growth Rate", SEmult = 2, bgcolor = 0) + # add bgcolor = 0 if you want a b/w background
  title(main = "Chr 5, Pos 76.6 cM") # Change manually so that pos & chr match above


#####################################################################################################
#################################### PERCENT VARIANCE EXPLAINED #####################################
#####################################################################################################
## Goal: combine all "g" lists into a dataframe (with labelled column name) for each marker
## Have another dataframe with phenotypes (CT$pheno)
## Create separate dataframes from g combo so that RILs w/ LL & that pheno are in separate columns from RILs with CC & pheno
### Basically want to do what plot_pxg does, but in dataframes instead of plots

# More thoughts:
## Create list of chr & pos values you can pull from to make a for loop
## for loop will use maxmarg() to make a dataframe for allele types at each marker

# Notes:
## g = list of alleles by RIL at specific marker, need to give chr & pos
## geno = dataframe of chr and position locations of my peaks for each phenotype
## pheno = dataframe of all phenotypes by RIL from CT$pheno (matrix)
library(tidyverse)
library(dplyr)
peakLocations = read.csv("~/Desktop/QTL_Maps/PeaksData/2021-05-09_peaks.csv")

g <- maxmarg(pr, map, chr = 1, pos = 21.509791, return_char = TRUE, minprob = 0)

pheno = data.frame(CT$pheno)

# For loop to combine all allele data for each RIL at the peak marker location
geno = data.frame(1:162) # create empty dataframe
for (i in 1:nrow(peakLocations)) {
  print(peakLocations$lodcolumn[i])
  print(peakLocations$chr[i])
  print(peakLocations$pos_cM[i])
  g <- maxmarg(pr, map, chr = peakLocations$chr[i], pos = as.double(peakLocations$pos_cM[i]), return_char = TRUE, minprob = 0)
  geno = cbind(geno, g)
}

# Saved this work
## Manually made column names with marker chr & pos b/c couldn't figure that out here
# write.csv(geno, "~/Desktop/QTL_Maps/PeaksData/2021-05-09_allelesAtAllPeaks.csv")

# Re-read in the geno file (with manually fixed column names!)
geno = read_csv("~/Desktop/QTL_Maps/PeaksData/2021-05-09_allelesAtAllPeaks.csv", rownames_included = FALSE)

# Make dataframe with marker allele data & phenotype, by RIL
## Starting with growthRate at Chr 1, Pos 21.5 peak
data = cbind(geno, pheno)
peaks = colnames(data)
##########################################################################################
# The following is for one peak at a time only [5-14-2021]
## Need a for loop (or function) to do this for the others

# Workshopping for loop
# column = "column name"
# trait = "trait name"

# Gather posVmax data for the RILs with CC allele at chosen marker
growthRate_CCdata  = data %>%
  select(RIL, chr1_pos21.50979133_growthRate, overallGrowthRate) %>%
  filter(chr1_pos21.50979133_growthRate == "CC")
# Gather posVmax data for the RILs with LL allele at chosen marker
growthRate_LLdata  = data %>%
  select(RIL, chr1_pos21.50979133_growthRate, overallGrowthRate) %>%
  filter(chr1_pos21.50979133_growthRate == "LL")

# Calc means & var for posVmax for each allele
avgCC = mean(growthRate_CCdata$overallGrowthRate, na.rm = TRUE) # avg = 665.0547
avgLL = mean(growthRate_LLdata$overallGrowthRate, na.rm = TRUE) # avg = 594.3407
avgDiff = avgCC-avgLL # Think about whether you need absolute value here...
# avgDiff = avgLL-avgCC # Use if LL data is higher than CC data!
print(avgDiff)
var_overallGrowthRate = var(data$overallGrowthRate, y = NULL, na.rm = TRUE) # var = 6284.778


# Use avgDiff to subtract that much from higher allele to get same population group
adj_growthRate_LLdata = growthRate_LLdata %>%
  mutate(adj_growthRate = overallGrowthRate + avgDiff) # Plus or minus? PAY ATTENTION!!! 

# Now need to add adjusted Vmax of LL to CC & make new posVmax column to calc var from
# adjusted_overallGrowthRate = combine(adj_overallGrowthRate_LLdata$adj_overallGrowthRate, overallGrowthRate_CCdata$overallGrowthRate)
adjusted_overallGrowthRate = append(adj_growthRate_LLdata$adj_growthRate, growthRate_CCdata$overallGrowthRate)



adj_var_overallGrowthRate = var(adjusted_overallGrowthRate, y = NULL, na.rm = TRUE) # var = 5057.851
percentVarianceExplained = 1-(adj_var_overallGrowthRate/var_overallGrowthRate) # perVarExplained = 0.1952219
##########################################################################################


# Max & min of posVmax
max(data$positionVmax, na.rm = TRUE) # max = 846.875
min(data$positionVmax, na.rm = TRUE) # min = 384.4

# Calc overall variance for each phenotype
var_Vmax = var(pheno$Vmax, y = NULL, na.rm = TRUE) # var = 1.698002e-07
var_vf = var(pheno$vf, y = NULL, na.rm = TRUE) # var = 0.3065314
var_positionVmax = var(pheno$positionVmax, y = NULL, na.rm = TRUE) # var = 6284.778
var_growthZoneWidth = var(pheno$growthZoneWidth, y = NULL, na.rm = TRUE) # var = 20661.62



#####################################################################################################
##### Allele effect (Estimated QTL Effect) #####
# scan1coef() - This function takes a single phenotype and the genotype probabilities 
##              for a single chromosome and returns a matrix with the estimated  
##              coefficients at each putative QTL location along the chromosome.
# (pr = 1) = chr 1
Aeff <- scan1coef(pr[,"4"], CT$pheno[,"growthZoneWidth"])
par(mar = c(4.1, 4.1, 1.1, 2.6), las = 1)
col <- c("slateblue", "violetred", "green3")
# map["2"] means chr 2, columns = 1:3 b/c 3 columns in Aeff (L, C, and intercept)
plot(Aeff, map["4"], columns = 1:3, col = col)
last_coef <- unclass(Aeff)[nrow(Aeff),] # pull out last coefficients
for(i in seq(along = last_coef))
    axis(side = 4, at = last_coef[i], names(last_coef)[i], tick = FALSE, col.axis = col[i])


##########################################################################################
# SNP effect or QTL effect, chromosome 5 and trait T50
c5eff <- scan1coef(pr[,"5"],CT$pheno[,"T50"]) 
#par(mfrow=c(2,1), mar=c(5,2,2,2), pin=c(3,3))
#plot_coef(c5eff[,1:2],map)
plot_coef(c5eff[,1:2],map,ylim=c(-1.5, 1.5))
#or
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
col <- c("slateblue", "violetred")
plot(c5eff, map["5"], columns=1:2, col=col)
last_coef <- unclass(c5eff)[nrow(c5eff),] # pull out last coefficients
for(i in seq(along=last_coef))
    axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])






#mike broke this, its a plotting loop which goes through your phenotypes and populations. fix later 
 for (i in 1:ncol(out1)){
    pdf(paste0("out/",colnames(out1)[i],'.pdf'), width=11, height=3)
    plot(out1, CT$gmap, lodcolumn = i, ylim=c(0,maxlod(out1*1.05)), main=paste0(colnames(out1)[i]))
    abline(h=sig, col="red")
    dev.off()
  }  

  pdf(paste0(setname,'.pdf'), width=11, height=3)
  plot(out1, CT$gmap, ylim=c(0,7),main=setname)
  abline(h=sig, col="red")
  dev.off()
  