\name{mam_sparse_dr}
\alias{mam_sparse_dr-function}
\alias{mam_sparse_dr}
\docType{package}
\title{
  Fit MAM with sparsity assumption and ranks selected by \code{BIC} or \code{CV}.
}
\description{
  Fit a high-dimensional multivariate additive model using B-splines, with or with aparsity assumptions and ranks selected by \code{BIC} or \code{CV}. The multivariate sparse group \code{lasso} (\code{mcp} or \code{scad}) and the coordinate descent algorithm are used to estimate
  functions for sparsity situation. The tuning parameter is selected by \code{BIC} or \code{CV}, which matchs the method of rank selection.
}

\usage{
mam_sparse_dr(Y, X, method = "BIC", ncv = 10, penalty = "LASSO", K_index = NULL,
              r1_index = NULL, r2_index = NULL, r3_index = NULL, lambda = NULL,
              SABC = NULL, nlam = 50, degr = 3, lam_min = 0.01,
              eps1 = 1e-4, maxstep1 = 20, eps2 = 1e-4, maxstep2 = 20, gamma = 2,
              dfmax = NULL, alpha = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{A \code{\eqn{n\times q}} numeric matrix of responses.}
  
  \item{X}{A \code{\eqn{n\times q}} numeric design matrix for the model.}
  
  \item{method}{The method to be applied to select parameters.  Either "BIC"
                (the default), or "CV".}
    
  \item{ncv}{The number of cross-validation folds.  Default is 10. If \code{method} is \code{BIC}, \code{ncv} is useless.}
  
  \item{penalty}{The penalty to be applied to the model. Either "LASSO" (the default), 
                 "SCAD", or "MCP".}
                 
  \item{K_index}{A user-specified sequence of \code{K} values, where \code{K} is he number 
                 of B-spline base function. Default is \code{k_index=6}.}
                 
  \item{r1_index}{A user-specified sequence of \code{\eqn{r_1}} values, where 
                  \code{\eqn{r_1}} is the first dimension of single value matrix of the  
                  tensor. Default is 
                  \code{r1_index\eqn{=1,\cdots,\min(\lceil\log(n)\rceil,p)}}.}
            
  \item{r2_index}{A user-specified sequence of \code{\eqn{r_2}} values, where 
                  \code{\eqn{r_2}}is the second dimension of single value matrix of the
                  tensor. Default is\code{r2_index\eqn{=1,\cdots},max{K_index}.}
                  }
                  
  \item{r3_index}{A user-specified sequence of \code{\eqn{r_3}} values, where 
                  \code{\eqn{r_3}} is the third dimension of single value matrix of the  
                  tensor. Default is 
                  \code{r3_index\eqn{=1,\cdots,\min(\lceil\log(n)\rceil,q)}}.}
                  
  \item{lambda}{A user-specified sequence of lambda values.  By default,
                a sequence of values of length \code{nlam} is computed, equally
                spaced on the log scale.}

  \item{SABC}{A user-specified list of initial coefficient matrix of \code{\eqn{S}},  
              \code{\eqn{A}}, \code{\eqn{B}}, \code{\eqn{C}}. By default,
              initial matrices are provided by random.}
              
  \item{nlam}{The number of lambda values.  Default is 50.}
  
  \item{degr}{The number of knots of B-spline base function. Default is \code{degr = 3}.}
  
  \item{lam_min}{The smallest value for lambda, as a fraction of
                 lambda.max.  Default is 1e-2.}
                 
  \item{eps1}{Convergence threshhold.  The algorithm iterates until the
              relative change in any coefficient is less than \code{eps1}.  
              Default is \code{1e-4}.}
              
  \item{maxstep1}{Maximum number of iterations. Default is 20.}
  
  \item{eps2}{Convergence threshhold.  The Coordinate descent method algorithm iterates 
              until the relative change in any coefficient is less than \code{eps2}.  
              Default is \code{1e-4}.}
  
  \item{maxstep2}{The maximum iterates number of coordinate descent method.  Default is 20.}
  
  \item{gamma}{The tuning parameter of the MCP/SCAD penalty (see details).}
  
  \item{dfmax}{Upper bound for the number of nonzero coefficients.
               Default is no upper bound.  However, for large data sets,
               computational burden may be heavy for models with a large number of
               nonzero coefficients.}
    
  \item{alpha}{Tuning parameter for the Mnet estimator which controls
               the relative contributions from the LASSO, MCP/SCAD penalty and the ridge,
               or L2 penalty.  \code{alpha=1} is equivalent to LASSO, MCP/SCAD penalty,
               while \code{alpha=0} would be equivalent to ridge regression.
               However, \code{alpha=0} is not supported; \code{alpha} may be
               arbitrarily small, but not exactly 0.}

}

\details{
  This function gives \code{pq} functional coefficients' estimators of MAM. The singular value matrix of 
  tensor is a \code{\eqn{r_1\times r_2\times r_3}}-tensor. We choose \code{\eqn{r_1}}, \code{\eqn{r_2}} 
  and \code{\eqn{r_3}}  by \code{BIC} or \code{CV}.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use

  \item{Dnew}{Estimator of \code{\eqn{D_{(3)}}}.}
  
  \item{rss}{Residual sum of squares (RSS).}
  
  \item{df}{Degrees of freedom.}
  
  \item{activeA}{The active set.}
  
  \item{lambda}{The sequence of regularization parameter values in the path.}
  
  \item{selectedID }{The index of \code{lambda} corresponding to \code{lambda_opt}.}
  
  \item{lambda_opt }{The value of \code{lambda} with the minimum \code{BIC} or \code{CV} value.}
  
  \item{RSS}{The values of \code{BIC} or \code{CV}, which is a vector.}
  
  \item{rk_opt}{The optimal parametres that slected by \code{BIC} or \code{CV}. It is a vector with length 4, which are selected \code{\eqn{r_1},  \eqn{r_2}, \eqn{r_3}}, and \code{\eqn{K}}.}
  
  \item{Y}{Response \code{\eqn{Y}}.}
  
  \item{X}{Design matrix \code{\eqn{X}}.}
  
  \item{Z}{Design matrix of Bspline approximation.}
  %\item{...}{ Other options for CompositeQuantile.}
}

%\author{
%Your Name, email optional.
%Maintainer: Xu Liu <liu.xu@sufe.edu.cn>
%}
\references{
  A tensor estimation approach to multivariate additive models.
}
\seealso{
  mam_dr, mam_sparse
}

\examples{
  #Example 1

  D2 <- matrix(runif(50, 0.7, 1), 2, 25) # tensor with size 5*2*5
  mydata <- generateData(200, 5, 10, 5, D2)
  
  fit <- mam_sparse_dr(mydata$Y, mydata$X)
  D3hat <- fit$Dnew
  D2hat <- TransferModalUnfoldings(D3hat,3,2,5,2,5)
  opt <- fit$rk_opt

  #Example 2

  data(breastdata.rda)
  Y = t(breastdata$dna[breastdata$chrom==21,]) 
  Xt = t(breastdata$rna[which(breastdata$genechr==21),])
  n = nrow(Y)

  minX = apply(Xt,2,min)
  maxX = apply(Xt,2,max)
  X = (Xt - matrix(rep(minX,each = n),n))/matrix(rep(maxX-minX,each = n),n)
  Y = scale(Y)
  fit <- mam_sparse_dr(Y, X)
  D3hat <- fit$Dnew
  opt <- fit$rk_opt
}