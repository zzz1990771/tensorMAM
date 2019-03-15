# tensorMam
 A tensor estimation approach to multivariate additive models.
 
  For a high-dimensional multivariate additive model (MAM) using B-splines, with or without aparsity assumptions, 
  treating the coefficients as a third-order tensor and borrowing Tucker decomposition to reduce the number of parameters.  
  The multivariate sparse group lasso (mcp or scad) and the coordinate descent algorithm are used to estimate
  functions for sparsity situation.
# Installation

    #install.packages("devtools")
    #install Rtools 3.5 (http://cran.r-project.org/bin/windows/Rtools)
    library(devtools)
    install_github("xliusufe/tensorMam")

# Usage

   - [x] [tensorMam-manual](https://github.com/xliusufe/tensorMam/blob/master/inst/tensorMam-manual.pdf) ------------ Details of the usage of the package.
# Example

    library(tensorMam)

    D2 <- matrix(runif(50, 0.7, 1), 2, 25)
    mydata <- generateData(200, 5, 5, 5, D2)    
    fit <- mam(mydata$Y, mydata$X)
    coeff <- fit$Dnew
    
    fit_dr <- mam_dr(mydata$Y, mydata$X)
    opt <- fit_dr$rk_opt
 
 # References
A tensor estimation approach to multivariate additive models. Manuscript.

# Development
The R-package is developed by Xu Liu (liu.xu@sufe.edu.cn) and Xiangyong Tan.
