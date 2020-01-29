# tensorMam
 A tensor estimation approach to multivariate additive models.
 
  For a high-dimensional multivariate additive model (MAM) using B-splines, with or without aparsity assumptions, 
  treating the coefficients as a third-order tensor and borrowing Tucker decomposition to reduce the number of parameters.  
  The multivariate sparse group lasso (mcp or scad) and the coordinate descent algorithm are used to estimate
  functions for sparsity situation.
# Installation

    #install Rtools 3.5 (http://cran.r-project.org/bin/windows/Rtools)
    #install.packages("devtools")
    #install.packages("Rcpp")
    library(devtools)
    install_github("xliusufe/tensorMam")

# Usage

   - [x] [tensorMam-manual](https://github.com/xliusufe/tensorMam/blob/master/inst/tensorMam-manual.pdf) ------------ Details of the usage of the package.
# Example
    
    library(tensorMam)
    # Example 1
    # The usage of function "mam()"
    p <- 5
    q <- 5
    D2 <- matrix(runif(2*p*q, 0.7, 1), 2, p*q) 
    mydata <- generateData(200, q, p, p, D2) 
    fit <- mam(mydata$Y, mydata$X)
    K <- fit$K
    D3hat <- fit$Dnew
    D2hat <- TransferModalUnfoldings(D3hat,3,2,p,K,q)
    mu <- fit$mu
    
    # Example 2
    # The usage of function "mam_dr()"
    fit_dr <- mam_dr(mydata$Y, mydata$X)
    K <- fit_dr$K
    D3hat <- fit$Dnew
    D2hat <- TransferModalUnfoldings(D3hat,3,2,p,K,q)	
    mu <- fit$mu
    opt <- fit_dr$rk_opt	
    
    # Example 3 
    # The usage of function "mvrblockwise()"
    n <- 200
    q <- 5
    s <- 3
    p <- 100
    B <- matrix(runif(q*s, 2,3), s)
    X <- matrix(rnorm(n*p),n,p)
    Y <- X[,1:s]%*%B + matrix(rnorm(n*q),n)
    fit <- mvrblockwise(Y,X) #See details in the function "mvrblockwise"
    fit$activeX
    fit$Bhat
    which(rowSums(fit$Bhat^2)>0)
    fit$muhat
    
    # Example 4
    # The usage of function "mvrcolwise()"
    n <- 200
    q <- 5
    s <- 3
    p <- 100
    B <- matrix(runif(q*s, 2,3), s)
    X <- matrix(rnorm(n*p),n,p)
    Y <- X[,1:s]%*%B + matrix(rnorm(n*q),n)
    fit <- mvrcolwise(Y,X) #See details in the function "mvrcolwise"
    fit$activeX
    fit$Bhat
    which(rowSums(fit$Bhat^2)>0)
    fit$muhat
    
 
 # References
A tensor estimation approach to multivariate additive models. Manuscript.

# Development
The R-package is developed by Xu Liu (liu.xu@sufe.edu.cn) and Xiangyong Tan.
