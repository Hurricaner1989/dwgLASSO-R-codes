# A demo of dwgLASSO to prioritize gene list for surival time prediction 
#
# Reference:
# [1] Zuo, Yiming, Yi Cui, Guoqiang Yu, Ruijiang Li, and Habtom W. Ressom. "Incorporating prior 
#     biological knowledge for network-based differential gene expression analysis using 
#     differentially weighted graphical LASSO." BMC Bioinformatics 18, no. 1 (2017): 99.
#
# Copyright 2015-2017, Yiming Zuo.

## Start clean
rm(list = ls())

## Make sure current working directory is dwgLASSO package
getwd()

## Load data: one microarray dataset with the survival records of patients
# The dataset has been properly preprocessed
# X represents the expression value of 58 genes on 158 patients
# y represents the survival time (months) of the 158 patients
# c represents whether the patient is alive (0) or dead (1) up to the observation time
load("Data/Data.RData")

## Load R packages
# glasso 
library("glasso") # need install 
# mvtnorm to compute log likelihood error
library("mvtnorm") # need install

## Load functions
source("function.R")

## Sort genes based on their adjusted p-values
pvalue_adjust_idx <- sort(pvalue.adjust, index.return = T)$ix
sig_l_s <- data.frame("gene symbol" = gene.list[pvalue_adjust_idx],
                      "adjusted pvalue" = pvalue.adjust[pvalue_adjust_idx],
                      stringsAsFactors = F)
rm(gene.list) # delete useless variable

## Create W matrix 
# save all significant genes to txt
write.table(sig_l_s$gene.symbol, file = "Data/sig_gene.txt", quote = F, row.names = F, 
            col.names = F)
#################################################################################################
# input the list to String database (http://string-db.org/) and get the confidence score saved as
# sig_string.txt in Data folder
#################################################################################################
# read in the confidence score
w_list <- read.table("Data/sig_string.txt", header = TRUE, sep = "\t", stringsAsFactors = F)
# initialize weight matrix
W <- matrix(0, nrow(sig_l_s), nrow(sig_l_s))
# test the range of the combined score, min is 0.15 in this list
#max(w_list$combined_score) 
#min(w_list$combined_score)
# assign confidence score to W
for (k in 1:nrow(w_list)) {
    i <- match(w_list[k,]$node1, sig_l_s$gene.symbol)
    j <- match(w_list[k,]$node2, sig_l_s$gene.symbol)
    W[i,j] <- w_list[k,]$combined_score
    W[j,i] <- W[i,j]
}
sum(W!=0) # test how many entries have non-zero value, 238 in this list
rm(i, j, k, w_list)

## Divide data into low risk and high risk groups
# define low risk (data_low: survival time more than 60 months) and high risk groups 
# (data_high: survival time less than 60 months and dead)
data <- group(pvalue_adjust_idx, X, c, y)
rm(X, y, c)

## Z-transform the data for group-specific normalization
data_low <- scale(t(data$low.risk.group))
data_high <- scale(t(data$high.risk.group))
cov_low <- var(data_low)
cov_high <- var(data_high)
rm(data)

set.seed(1)
## Apply weighted graphical LASSO (wgLASSO) given W, data_z_low and data_z_high
# low risk first
n_fold <- 10 # number of folds
rho = exp(seq(log(0.4), log(0.01), length.out = 20))
Ones <- matrix(rep(1, ncol(data_low)), nrow = ncol(data_low), ncol = ncol(data_low))
# draw error curve 
error <- choose_rho(data_low, n_fold, rho)
# chosse optimal rho 
rho[error$log.cv == min(error$log.cv)] # rho based on minimum rule
abline(v = rho[error$log.cv == min(error$log.cv)], col = "red", lty = 3)
# one standard error rule
abline(h = min(error$log.cv) + error$log.rho[error$log.cv == min(error$log.cv)], col = "blue") 
rho
#################################################################################################
# select rho that are under blue line and on the right of the red line, in this case, rho = 0.184
#################################################################################################
rho_opt <- 0.184
# perform wgLASSO
pre_low <- glasso(cov_low, rho = as.matrix(rho_opt * (Ones - W)))
rm(n_fold, rho, rho_opt, Ones, error)
# covert numeric prcision matrix to binary matrix
thres <- 1e-3
pre_low_binary <- pre_low$wi
pre_low_binary[abs(pre_low$wi) < thres] = 0
pre_low_binary[abs(pre_low$wi) > thres] = 1
diag(pre_low_binary) = 0
sum(pre_low_binary != 0) # how many connections in network, 768
rm(thres)

# high risk second
n_fold <- 10 # number of folds
rho = exp(seq(log(0.4), log(0.01), length.out = 20))
Ones <- matrix(rep(1, ncol(data_low)), nrow = ncol(data_low), ncol = ncol(data_low))
# draw error curve 
error <- choose_rho(data_high, n_fold, rho)
# chosse optimal rho 
rho[error$log.cv == min(error$log.cv)] # rho based on minimum rule
abline(v = rho[error$log.cv == min(error$log.cv)], col = "red", lty = 3)
# one standard error rule
abline(h = min(error$log.cv) + error$log.rho[error$log.cv == min(error$log.cv)], col = "blue") 
rho
#################################################################################################
# select rho that are under blue line and on the right of the red line, in this case, rho = 0.223
#################################################################################################
rho_opt <- 0.223
# perform wgLASSO
pre_high <- glasso(cov_high, rho = as.matrix(rho_opt * (Ones - W)))
rm(n_fold, rho, rho_opt, Ones, error)
# covert numeric prcision matrix to binary matrix
thres <- 1e-3
pre_high_binary <- pre_high$wi
pre_high_binary[abs(pre_high$wi) < thres] = 0
pre_high_binary[abs(pre_high$wi) > thres] = 1
diag(pre_high_binary) = 0
sum(pre_high_binary != 0) # how many connections in network, 760
rm(thres)

## Network analysis: dwgLASSO
# node degree
degree_low <- rowSums(pre_low_binary) 
degree_high <- rowSums(pre_high_binary)
# scaled node degree
degree_low_s <- (degree_low - min(degree_low)) / (max(degree_low) - min(degree_low)) 
degree_high_s <- (degree_high - min(degree_high)) / (max(degree_high) - min(degree_high))
# rank based on differential network score
dns <- abs(degree_low_s - degree_high_s)
degree_diff_idx <- sort(dns, decreasing = T, index.return=T)$ix
sig_l_s$gene.symbol[degree_diff_idx[1:10]] # top 10 genes based on dwgLASSO
sig_l_s$adjusted.pvalue[degree_diff_idx[1:10]] # their adjusted p-values 

## save top 10 significant genes after prioritization based on dwgLASSO
write.table(sig_l_s[degree_diff_idx[1:10],], file = "sig10table_dwgLASSO.csv", sep=",", quote = F, 
            row.names = F, col.names = F)
