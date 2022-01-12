QTLCode.R
library(data.table)
library(qtl2)
useQTLCode = function(){
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
  CT <- read_cross2("~/Desktop/QTL_AnalysisMaterials/Ashley_QTL/CvixLer.workingTMP.json")
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
  out <- scan1(pr, CT$pheno[,c(9)], k_loco, cores = 4) ## Note: (1:4) for REGR Curve Descriptors, (1,2,8,9) for parameters
  
  
  summary(out)
  nperm = 2
  # Permutation testing
  operm <- scan1perm(pr[,1:5], CT$pheno[,c(9)], k_loco[1:5], addcovar = NULL,
                     Xcovar = NULL, intcovar = NULL, weights = NULL, reml = TRUE,
                     model = "normal", n_perm = nperm, perm_Xsp = FALSE, # changed n_perm=10000 from 1000 (5.9.2021)
                     perm_strata = NULL, chr_lengths = NULL, cores = 4)
  
  # sig <- summary(operm, alpha = c(0.05, 0.01)) # Significant threshold, 95% & 99%  
  sig <- summary(operm, alpha = c(0.05)) # Significant threshold, just 95%
  # ymx <- maxlod(out) # max lod score
  
  ##### Finding LOD peaks
  # 
  peaks <- find_peaks(out, map, threshold = sig[1,1], drop = 1.5)
  # peaks_fullDataset <- find_peaks(out, map, threshold = sig[1,1], drop = 1.5)
  # peaks_noMatZone <- find_peaks(out, map, threshold = sig[1,1], drop = 1.5)
  # peaks_littleMatZone <- find_peaks(out, map, threshold = sig[1,1], drop = 1.5)
  
  peaks
  
}