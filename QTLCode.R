# install.packages("data.table")
# install.packages("qtl2", repos="http://rqtl.org/qtl2cran")

library(data.table)
library(qtl2)

# Folder with all of the input files for this R script
# setwd("~/ashleyhenry/Desktop/QTL_Wu/SamebutwithCviLer")
setwd("~/Desktop/QTL_WuPracticewithSteveDeslauriers/Ashley_QTL")

##################################################################################################
########################## Use this chunk if need to add leading zeros ###########################
############## Tried on 2020-05-27 & royally ruined everything for the QTLCode.R! ################
##################################################################################################
# Need to add leading zeros to genotype data: CvixLer.RIL.numSNPs.csv
## Load in data, add zeros to be 3 digits after RIL
## Save new RIL numbers with leading zeros as csv file in working directory to use in QTL map making
# geno <- read.csv("~/Desktop/QTL_WuPracticewithSteveDeslauriers/CvixLer/Ashley_QTL/CvixLer.RIL.SNPs.csv")
# geno$Genotype <- sprintf("RIL%03d", geno$Genotype)
# geno %>%
#   order(geno$Genotype)
# write.csv(geno, "CvixLer.RIL.SNPs.csv")
# welp this now says it can't read either geno or pheno files... back to the olde drawing board...
##################################################################################################

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
  CT <- read_cross2("~/Desktop/QTL_WuPracticewithSteveDeslauriers/Ashley_QTL/CvixLer.working.json")
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
  k_loco <- calc_kinship(pr, "loco") # produced list of lists containing #s from 0.7-0.9 (what is this??)
# QTL mapping. cols: 1=x0, 2=vf, 3=k, 4=n (I think) 

## Use when selecting a range of columns, make sure to follow through in scan1perm() below
  # out <- scan1(pr, CT$pheno[,2:6], k_loco, cores = 4) 
# Use when selecting specific columns
  out <- scan1(pr, CT$pheno[,c(1,2,5,6)], k_loco, cores = 4) 
  # Attempting scantwo function, as of 9.9 it's not going well...
  # out2 <- scantwo(pr, CT$pheno[,c(1:4)], k_loco, verbose=FALSE)
  # 
  # plot(out2, chr=c(1:5))
  summary(out)
  
# Permutation testing
  operm <- scan1perm(pr[,1:5], CT$pheno[,c(1,2,5,6)], k_loco[1:5], addcovar = NULL,
                   Xcovar = NULL, intcovar = NULL, weights = NULL, reml = TRUE,
                   model = "normal", n_perm = 1000, perm_Xsp = FALSE,
                   perm_strata = NULL, chr_lengths = NULL, cores = 4)

  sig <- summary(operm, alpha = c(0.05, 0.01)) # Significant threshold               
# ymx <- maxlod(out) # max lod score

##### Finding LOD peaks
  peaks <- find_peaks(out, map, threshold = sig[1,1], drop = 1.5)
  peaks

##### QTL Mapping Plots: lodcolumn 1 = x0; lodcolumn 2 = vf; lodcolumn 3 = k; lodcolumn 4 = n
  plot(out, map, lodcolumn = 1, ylim = c(0,8), col = rgb(0, 0, 0, 0.5), main = "QTL Map for Average Parameter Data") # vf = black
  plot(out, map, lodcolumn = 2, col = rgb(1, 0, 0, 0.5), add = TRUE) # log(k) = red
  plot(out, map, lodcolumn = 3, col = rgb(0, 0, 1, 0.5), add = TRUE) # log(n) = blue
  plot(out, map, lodcolumn = 4, col = rgb(0, 1, 0, 0.5), add = TRUE) # n = green
  # plot(out, map, lodcolumn = 5, col = rgb(0.4, 0, 0.6, 0.5), add = TRUE)# log(n) = violet

  abline(h = sig[,1], col = rgb(0, 0, 0, 0.5)) # vf = black
  abline(h = sig[,2], col = rgb(1, 0, 0, 0.5)) # log(k) = red
  abline(h = sig[,3], col = rgb(0, 0, 1, 0.5)) # log(n) = blue
  abline(h = sig[,4], col = rgb(0, 1, 0, 0.5)) # n = green
  # abline(h = sig[,5], col = rgb(0.4, 0, 0.6, 0.5)) # log(n) = violet
  legend("topright", lwd = 3, col = c(rgb(0, 0, 0, 0.5), rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5), rgb(0, 1, 0, 0.5), rgb(0.4, 0, 0.6, 0.5)), colnames(out), bty = "n")

plot(out, map[peaks[2,3]], lodcolumn = 2, xlim = c(peaks[2,6],peaks[2,7]),
                                      ylim = c(0,7), main = paste0(colnames(out)[2]))
dev.off()

#####################################################################################################
###### Added by Ashley: Checking the marker regression at markers with the highest LOD scores #######
#####################################################################################################
# Plots geno v pheno to see data spread 
## Compares averages between parent genotypes at that cM position
par(mfrow = c(2,2))
# Data = out, markers = ""
plotPXG(out, "")
plotPXG(out, "")

#####################################################################################################
##### Allele effect #####
Aeff <- scan1coef(pr[,"5"], CT$pheno[,"T50"])
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
col <- c("slateblue", "violetred", "green3")
plot(Aeff, map["5"], columns=1:3, col=col)
last_coef <- unclass(Aeff)[nrow(Aeff),] # pull out last coefficients
for(i in seq(along=last_coef))
    axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])


g <- maxmarg(pr, map, chr=5, pos=181255938, return_char=TRUE, minprob = 0)
#g <- maxmarg(pr, map, chr=10, pos=145648045, return_char=TRUE, minprob = 0)
#par(mar=c(4.1, 4.1, 1, 0.6))
plot_pxg(g, CT$pheno[,"T50cont"], ylab="T50 (days)", SEmult=1)
title(main="Allele effect at 181255938, chromosome 5")

####### Emergency rate #######
out1 <- scan1(pr, CT$pheno[,5], k_loco, cores=4) 
summary(out1)
operm1 <- scan1perm(pr[,1:10], CT$pheno[,5], k_loco[1:10], addcovar = NULL,
                   Xcovar = NULL, intcovar = NULL, weights = NULL, reml = TRUE,
                   model = "normal", n_perm = 1000, perm_Xsp = FALSE,
                   perm_strata = NULL, chr_lengths = NULL, cores = 4) 
sig1 <- summary(operm1, alpha=c(0.05, 0.01)) # Significant threshold  

peaks1 <- find_peaks(out1, map, threshold=sig1, drop=1.5)
peaks1

plot(out1,map,lodcolumn=1,ylim=c(0,7),col="slateblue", main="Seedling Emergence Rate")
abline(h=sig1, col="violetred")

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
