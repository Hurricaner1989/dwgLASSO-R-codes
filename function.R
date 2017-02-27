# Functions defined to help demo.R
#
# Reference:
# [1] Zuo, Yiming, Yi Cui, Guoqiang Yu, Ruijiang Li, and Habtom W. Ressom. "Incorporating prior 
#     biological knowledge for network-based differential gene expression analysis using 
#     differentially weighted graphical LASSO." BMC Bioinformatics 18, no. 1 (2017): 99.
#
# Copyright 2015-2017, Yiming Zuo.

## Create group function to divide data into low risk and high risk groups
group <- function(pvalue_idx, X, c, y) {
    X_pn <- t(X)[pvalue_idx,]
    data_low <- X_pn[, which(y > 60)]
    data_high <- X_pn[, which((c == 1) & (y < 60))]
    data_list <- list("low.risk.group" = data_low, "high.risk.group" = data_high)
    return(data_list)
}

## Create log likelihood error function 
loglik_ave <- function(data, theta){
    loglik <- c()
    loglik <- log(det(theta))-sum(diag(var(data)%*%theta))
    -loglik
}

## Draw error curve
choose_rho <- function(data, n_fold, rho) {
    # randomly shuffle the data
    Data_low <- data[sample(nrow(data)),]
    # create n_fold equally size folds
    folds <- cut(seq(1, nrow(Data_low)), breaks = n_fold, labels = FALSE)
    # tune parameters
    d <- ncol(Data_low)
    Ones <- matrix(rep(1,d), nrow = d, ncol = d)
    loglik_cv <- c()
    loglik_rho <- c()
    pb <- txtProgressBar(min = 0, max = length(rho), style = 3) # create progress bar
    for(i in 1:length(rho)){
        Sys.sleep(0.1)
        # perform n_fold cross validation
        loglik <- c()
        for(j in 1:n_fold){
            # segement your data by fold using the which() function 
            testIndexes <- which(folds == j, arr.ind=TRUE)
            testData <- Data_low[testIndexes, ]
            trainData <- Data_low[-testIndexes, ]
            # use test and train data partitions however you desire...
            mu <- apply(t(trainData), 1, mean)
            cov <- var(trainData) # compute the covariance matrix
            pre <- glasso(cov, rho = as.matrix(rho[i] * (Ones - W)))
            loglik <- c(loglik, loglik_ave(testData, mu, pre$w))  
        }
        loglik_cv <- c(loglik_cv, sum(loglik) / n_fold)
        loglik_rho <- c(loglik_rho, sd(loglik) / sqrt(n_fold))
        setTxtProgressBar(pb, i) # update progress bar
    }
    close(pb)
    plot(rho, loglik_cv, main = "Error curve using corss validation",
         xlab = expression(lambda), ylab = "Error")
    lines(rho, loglik_cv)
    error <- list("log.cv" = loglik_cv, "log.rho" = loglik_rho)
    return(error)
}

