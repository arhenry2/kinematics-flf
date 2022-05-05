
# This is a new QTLCode to use scantwo() to find interacting qtl
## Had to make this to work with R/qtl as scantwo isn't in R/qtl2 (as far as I know?)
# Author: Ashley Henry
# Date: 23 March 2022

# Trying out Interaction from Broman's FAQ page
## 7 March 2022
library(qtl)
data(hyper)
hyper <- calc.genoprob(hyper, step=2.5, err=0.001)

qtl <- makeqtl(hyper, chr=15, pos=18, what="prob")

out.i <- addqtl(hyper, qtl=qtl, formula=y~Q1*Q2, method="hk")
out.a <- addqtl(hyper, qtl=qtl, formula=y~Q1+Q2, method="hk")
plot(out.i - out.a)

out2 <- scantwo(hyper, verbose=FALSE)
plot(out2, chr=c(1,4,6,7,15)) # pg. 217 of Broman's R/qtl manual
plot(out2, chr=c(1,4,6,7,15), lower="cond-int") # cleans up the plot to highlight pairs of QTL
plot(out2, chr=1, lower="cond-int", upper="cond-add")


# Trying it out on my data:
# Load in geno & pheno data (had to make new csvs from Broman's manual -> https://rqtl.org/sampledata/)
# root_data <- read.cross("csvs", "~/Desktop/", "ALFP_geno_new.csv", "2022-02-07_gr-nongr-maxREGR-width_new.csv)
# Judging by the scantwo() output, I'm trying this one phenotype at a time
## Didn't make it better... I don't know what I'm looking at...
root_data <- read.cross("csvs", "~/Desktop/", "ALFP_geno_new.csv", "2022-03-31_nonGR.csv")
# qtl <- makeqtl(hyper, chr=3, pos=18, what="prob")
# 
# out.i <- addqtl(hyper, qtl=qtl, formula=y~Q1*Q2, method="hk")
# out.a <- addqtl(hyper, qtl=qtl, formula=y~Q1+Q2, method="hk")
# plot(out.i - out.a)

root_data2 <- calc.genoprob(root_data, step=2.5, err=0.001)
out2 <- scantwo(root_data2, chr = c(3,4), verbose=FALSE) 
plot(out2, chr=c(3,4)) # pg. 217 of Broman's R/qtl manual ## full model vs null qtl
plot(out2, chr=c(3,4), lower="cond-int") # cleans up the plot to highlight pairs of QTL ## full model vs single
plot(out2, chr=c(3,4), lower="cond-int", upper="cond-add") ## single model vs null




root_data <- read.cross("csvs", "~/Desktop/", "ALFP_geno_new.csv", "2022-03-31_kineGR.csv")
root_data$pheno[,2] = CT$pheno[,c(6)]
root_data1 <- calc.genoprob(root_data, step=1, err=0.001)


out1 <- scanone(root_data1, chr = c(1:5), pheno.col = 2, verbose=FALSE)

# out1 <- scanone(root_data1, CT$pheno[,c(6)])
plot(out1, chr=c(1:5))





