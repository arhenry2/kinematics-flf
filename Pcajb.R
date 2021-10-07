# -*- coding: utf-8 -*-
# vim:fenc=utf-8
## ---------------------------
## Pcajb.R
## Purpose of script:
##
## Author: Julian Bustamante
##
## Date Created: 2021-08-04
##
## Copyright (c) Julian Bustamante, 2020
## Email: <jbustamante35@gmail.com>
##
## Distributed under terms of the MIT license.
##
## ---------------------------
##
## Notes:
## PCA class used for Ashley's tip angle data
##
## ---------------------------
library('pracma')

Pcajb <- setRefClass(
    "Pcajb",

    # Slots
    fields = list(
        Data = "matrix",
        npc  = "numeric"
    ),

    # Methods
    methods = list(
        meanSubtract = function() {
            x  <- Data
            u  <- colMeans(x)
            M  <- x - u
            Mu <- list(M,u)

            # Output Mu = [Mean-Subtracted Data] , [Column Means]
            return(Mu)
        },

        covarMatrix = function() {
            mu <- meanSubtract()
            m  <- Reshape(unlist(mu[1]), size(Data,1), size(Data,2))
            c  <- t(m) %*% m / size(m,1)

            # Output [Covariance Matrix]
            return(c)
        },

        eigens = function(neigs=0) {
            if (neigs == 0) {
                 neigs = npc
            }

            c  <- covarMatrix()
            wv <- eigen(c)

            # Eigen Vectors
            w  <- Reshape(unlist(wv[2]), size(c,1), size(c,2))
            w  <- w[,1:neigs]

            # Eigen Values
            v  <- Reshape(unlist(wv[1]), 1, size(c,1))
            v  <- v[1:neigs]

            wv <- list(w,v)

            # Output wv = [Eigen Vectors] , [Eigen Values]
            return(wv)
        },

        PCAScores = function(ndims=0,neigs=0) {
             if (ndims == 0) {
                 ndims = 1:size(Data,1)
             }

             if (neigs == 0) {
                 neigs = npc
             }

            x  <- Data[ndims,]
            mu <- meanSubtract()
            wv <- eigens(neigs)

            u  <- unlist(mu[2])
            w  <- Reshape(unlist(wv[1]), size(Data,2), neigs)
            s  <- pcaProject(x,w,u,'sim2scr')

            return(s)
        },

        SimData = function(ndims=0,neigs=0) {
            # Simulate data via projection and back-projection
             if (ndims == 0) {
                 ndims = 1:size(Data,1)
             }

             if (neigs == 0) {
                 neigs = npc
             }

            mu <- meanSubtract()
            wv <- eigens(neigs)
            s  <- PCAScores(ndims,neigs)

            u  <- unlist(mu[2])
            w  <- Reshape(unlist(wv[1]), size(Data,2), neigs)
            Y  <- pcaProject(s,w,u,'scr2sim')
            return(Y)
        },

        VarExplained = function(pct=0,neigs=0) {
            if (pct == 0) {
                pca = 1.0
            }

            if (neigs == 0) {
                neigs = npc
            }

            sz = size(Data,2)
            wv <- eigens(sz)
            v  <- Reshape(unlist(wv[2]), 1, sz)
            V  <- cumsum(v / sum(v))
            T  <- length(V[V <= pct]) + 1
            VT <- list(V[1:neigs],T)

            # Output VT = [Variance] , [Optimal PCs]
            return(VT)
        },

        pcaProject = function(x,w,u,req) {
            if (req == 'scr2sim') {
                # Back-project PC space to raw space
                m <- x %*% t(w)
                return(m + c(u))

            } else if (req == 'sim2scr') {
                # Project raw space to PC space
                m <- x - c(u)
                return(m %*% w)

            } else {
                return("Invalid projection")
            }
        }
    )
)
