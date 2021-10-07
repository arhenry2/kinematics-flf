# -*- coding: utf-8 -*-
# vim:fenc=utf-8
## ---------------------------
## testpca.R
##
## Author: Julian Bustamante
##
## Date Created: 2021-08-04
##
## Copyright (c) Julian Bustamante, 2021
## Email: <jbustamante@wisc.com>
##
## Distributed under terms of the MIT license.
##
## ---------------------------
##
## Notes:
##      Test my custom PCA class on Ashley's Data
##
## ---------------------------

setwd("/Users/ashleyhenry/kinematics-flf/")
library('pracma')
source('Pcajb.R')

# Load Data
N = 101
X = magic(N)

# Set PCA parameters
n  = 3
px = Pcajb(Data=X, npc=n)
z  = px$SimData()

# Plot Data
plot.new()
image(X)
image(z)
