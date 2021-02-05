## Following along with Example from Ch. 4 of Karl Broman's 
## A Guide to QTL Mapping with R/qtl
## Notes taken 21 May 2020

# Load in package and dataset
library(qtl)
data(hyper)

##### Book Section 4.1: Marker Regression #####

# Plots geno v pheno to see data spread 
## Compares averages between parent genotypes at that cM position
par(mfrow = c(1,2))
plotPXG(hyper, "D4Mit214")
plotPXG(hyper, "D12Mit20")

# Creates dataframe with marker, chr, cM position, and LOD score
## scanone() fits single-QTL models 
## Calcs LOD from several geno v pheno plots like those in line 12-13
## method = "mr" means marker regression (ANOVA at each marker)
out.mr <- scanone(hyper, method= 'mr')

# Quick threshold check to find LOD scores > 3
summary(out.mr, threshold = 3)
# Plots LOD scores from chr 4 & 12
## par() resets plotting area to one graph in frame
par(mfrow = c(1,1)) 
plot(out.mr, chr = c(4, 12), ylab = "LOD score")

##### Book Section 4.2: Interval Mapping #####

# Calculate the conditional genotype probabilities with calc.genoprob()
## Where: 
###  step = density (in cM) of the grid on which prob will be calc'd
### error.prob() = allows prob to be calc'd assuming given rate of genotyping errors
hyper <- calc.genoprob(hyper, step = 1, error.prob = 0.001)
hyper <- jittermap(hyper)
# Interval mapping using scanone() function; method = EM algorithm
out.em <- scanone(hyper, method = "em")
plot(out.em, ylab = "LOD score")
# Plotting both EM & MR method LOD scores on same plat
plot(out.em, out.mr, chr=c(4, 12), col=c("blue", "red"), ylab="LOD score")


# Notes from Ch. 8: 2-Dimensional, 2-QTL scans
library(qtl)
data(hyper)
# calc.genoprob() = calc the conditional QTL genotype prob, given the available marker data
# Using step = 2.5 as a course grid to increase computational speed (use higher for QTLCode.R?)
hyper <- calc.genoprob(hyper, step = 2.5, err = 0.001)
operm <- scanone(hyper, n.perm = 1000, verbose = FALSE) # this took 1.5 minutes to complete
# see how big hyper dataset is compared to yours

# scantwo() = 2-D, 2-QTL scan & default analysis is max likelihood via EM algorithm
# verbose=FALSE to reduce tracing informationw
out2 <- scantwo(hyper, verbose=FALSE)
# Plot LOD scores from 2-D scan for selected chr
plot(out2, chr=c(1,4,6,7,15))












