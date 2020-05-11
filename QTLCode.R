# install.packages("data.table")
# install.packages("qtl2", repos="http://rqtl.org/qtl2cran")

library(data.table)
#library(help=data.table)
library(qtl2)

# Folder with all of the input files for this R script
setwd("/Users/ashleyhenry/Desktop/QTL_Wu/SamebutwithCviLer")

############### Remove Tri-allelics, formating hmp to R/qtl ####################
############### recode snps into major minor 1:3 ###############################

# load file from Tassel5
map<-fread("B97XB73/B73xB97.104K.impute.BiAll.hmp.txt",stringsAsFactors=FALSE)

# replace IUPAC to bi-allelic dossage, set heterozygous to NA
map[map=="A"] <- "AA"
map[map=="G"] <- "GG"
map[map=="C"] <- "CC"
map[map=="T"] <- "TT"
map[map=="N"] <- NA
map[map=="R"] <- NA
map[map=="Y"] <- NA
map[map=="S"] <- NA
map[map=="W"] <- NA
map[map=="K"] <- NA
map[map=="M"] <- NA
map[map=="-"] <- NA
map[map=="0"] <- NA
#View (map)

# SNP list of tri-allelic from Tassel 5, filtered on Allele 3 freq in Excel
list <- read.csv("B97XB73/PruneList.csv", header=FALSE)
setlist <- as.character(list$V1)
str(setlist)
newdata <- map[!map$`rs#` %in% setlist,] # delete SNPs in PruneList from map
dim(newdata)

id<-newdata[,1:11]   # identification columns
SNP <- newdata[,12:ncol(newdata)] # subset file contains SNPs only
dim(SNP)

# scrime package to recode to 1,2,3 (Major homo het Minor homo)
RecodeSNP<-recodeSNPs(SNP, first.ref=FALSE) 
dim(RecodeSNP)
table(RecodeSNP[,5])
str(RecodeSNP)
done<-cbind(id, RecodeSNP) # add id columns

str(done)
mapinfo <-cbind(done$`rs#`, done$chrom, done$pos)
colnames(mapinfo) <- c("marker","chr","pos")
colnames(RecodeSNP)

rownames(RecodeSNP) <- done$`rs#`
# RecodeSNP[1:5,1:5]

inverseSNP<- t(RecodeSNP)

write.csv(inverseSNP, file="B97xB73.RIL.numSNPs.csv", row.names=FALSE, quote=FALSE)
write.csv(mapinfo, file="B97xB73.RIL.Map.csv", row.names = FALSE, quote = FALSE)

############# count total number of "A,T,G,C" for each TAXA ##################
############# pick one for each RIL with least missing data ##################
dat <- read.table("B73xB97.SNP.txt", header=T)
new <- data.frame(dat[,1])

# Count number of A/T/G/C per row and subtract by total number of columns
for (row in 1 : nrow(dat)) {
new[row,2] <- sum(dat[row,] == 'G'|dat[row,] == 'C'|dat[row,] == 'A'|dat[row,] == 'T')
  }
View(new)
write.csv(new, file="SNPsCount.csv")

############# RQTL on B97xB73 NAM, Cold Tolerance Project, MW, GW, CA, DL ######################
#/B97xB73    
  #/out                         -output folder from analysis
  #/B97xB73.working.json       - JSON which calls the objects dependents for rQTl2
    # /B97xB73.RIL.MAP.csv     - MAP information, with SNP name, CHR, and Pos for each sitename
    # /B97xB73.RIL.numSNPs.csv - SNPs in numeric coded format for rqtl, 
    # /B73xB97.phe.csv         - Phenotype.
################################### Mapping #####################################################
# read json object which is calling the MAP, NUM, and Pheno within rQTLmeta
  CT<-read_cross2("CvixLer.json")
  CT   # summary; crosstype "riself" for 2-way RIL by selfing 
  # names(CT)
  # head(CT$geno)
  # head(CT$pheno)
  g1 <- do.call("cbind", CT$geno) # CT$geno is a list of 10 chromosomes, combine 10 chromosomes
# Pseudomarkers improve resolution, but need huge memory. step=0 does nothing
  map <- insert_pseudomarkers(CT$gmap, step=0, stepwidth="max")
  pr <- calc_genoprob(CT, map, error_prob = 0.002)
# For each individual at each position, find the genotype with the maximum marginal probability.
  m1 <- maxmarg(pr)
# Estimate the numbers of crossovers in each individual on each chromosome.
  xo1 <- count_xo(m1)
  k_loco <- calc_kinship(pr, "loco")
# QTL mapping. col6 T50; col7 T50Cont
  out <- scan1(pr, CT$pheno[,1], k_loco, cores = 4) 
  summary(out)
# Permutation testing
  operm <- scan1perm(pr[,1:5], CT$pheno[,1], k_loco[1:5], addcovar = NULL,
                   Xcovar = NULL, intcovar = NULL, weights = NULL, reml = TRUE,
                   model = "normal", n_perm = 1000, perm_Xsp = FALSE,
                   perm_strata = NULL, chr_lengths = NULL, cores = 4) 

sig <- summary(operm, alpha=c(0.05, 0.01)) # Significant threshold               
#ymx <- maxlod(out) # max lod score

##### Finding LOD peaks
peaks <- find_peaks(out, map, threshold = sig[1,1], drop = 1.5)
peaks

##### T50 lodcolumn 1; T50Cont lodcolumn 2 colnames(out)[2]
plot(out, map, lodcolumn = 1, ylim = c(0,8), col = rgb(0, 0,1, 0.5), main = paste0(colnames(out)[1]))
plot(out,map,lodcolumn=2,col=rgb(1,0,0,0.5), main=paste0(colnames(out)[2]),add=TRUE)
abline(h=sig[,1], col=rgb(0,0,1,0.5)) # For T50
abline(h=sig[,2], col=rgb(1,0,0,0.5)) # For T50Cont
legend("topleft", lwd=2, col=c(rgb(0,0,1,0.5), rgb(1,0,0,0.5)), colnames(out), bty="n")

plot(out,map[peaks[2,3]],lodcolumn=2,xlim=c(peaks[2,6],peaks[2,7]),
                                      ylim=c(0,7), main=paste0(colnames(out)[2]))
dev.off()
##########################################################################################################
##### Allele effect
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
