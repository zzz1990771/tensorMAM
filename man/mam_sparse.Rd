\name{mam_sparse}
\alias{mam_sparse-function}
\alias{mam_sparse}
\docType{package}
\title{
  Fit MAM with sparsity assumption and fixed ranks.
}
\description{
  Fit a high-dimensional multivariate additive model using B-splines, with or without aparsity assumptions, and given ranks given ranks \eqn{r_1, r_2, r_3}. The multivariate sparse group \code{lasso} (\code{mcp} or \code{scad}) and the coordinate descent algorithm are used to estimate functions for sparsity situation.
}

\usage{
mam_sparse(y, x, K = 6, r1 = NULL, r2 = NULL, r3 = NULL, penalty = "LASSO",
           lambda = NULL, SABC = NULL, degr = 3, nlam = 20, lam_min = 1e-3, 
           eps1 = 1e-4, maxstep1 = 20, eps2 = 1e-4, maxstep2 = 20, gamma = 2,
           dfmax = NULL, alpha = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{A \code{\eqn{n\times q}} numeric matrix of responses.}
  
  \item{x}{A \code{\eqn{n\times q}} numeric design matrix for the model.}
  
  \item{K}{The number of B-spline base function, that is the plus of both degrees of base  
          functioin and the number of knots. Default is \code{7}.}
          
  \item{degr}{The number of knots of B-spline base function. Default is \code{3}.}
  
  \item{r1}{The first dimension of single value matrix of the tensor. Default is 2.}
  
  \item{r2}{The second dimension of single value matrix of the tensor. Default is 2.}
  
  \item{r3}{The third dimension of single value matrix of the tensor. Default is 2.}
  
  \item{penalty}{The penalty to be applied to the model. Either "LASSO" (the default),  
                 "SCAD", or "MCP".}
        
  \item{lambda}{A user-specified sequence of lambda values.  By default,
        a sequence of values of length \code{nlam} is computed, equally
        spaced on the log scale.}
        
  \item{SABC}{A user-specified list of initial coefficient matrix of \code{\eqn{S}},
              \code{\eqn{A}}, \code{\eqn{B}}, \code{\eqn{C}}. By default,
              initial matrices are provided by random.}
              
  \item{nlam}{The number of lambda values. Default is 20.}
   
  \item{lam_min}{The smallest value for lambda, as a fraction of
                 lambda.max.  Default is 1e-3.}
                 
  \item{eps1}{Convergence threshhold.  The algorithm iterates until the
              relative change in any coefficient is less than \code{eps1}.  
              Default is \code{1e-4}.}
              
  \item{maxstep1}{Maximum number of iterations. Default is \code{20}.}
  
  \item{eps2}{Convergence threshhold.  The Coordinate descent method algorithm iterates 
              until the relative change in any coefficient is less than \code{eps2}. 
              Default is \code{1e-4}.}
  
  \item{maxstep2}{The maximum iterates number of coordinate descent method.  
                  Default is \code{20}.}
                  
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
  \item{betapath}{Solution path of \code{\eqn{\beta}}.}
  
  \item{rss}{Residual sum of squares (RSS).}
  
  \item{df}{Degrees of freedom.}
  
  \item{lambda}{The sequence of regularization parameter values in the path.}
  
  \item{lambda_opt}{The value of \code{lambda} with the minimum
                    \code{BIC} value.}
  
  \item{selectedID}{The index of \code{lambda} corresponding to
                    \code{lambda_opt}.}
  
  \item{activeA}{The active set.}
  
  \item{Dnew}{Estimator of \code{\eqn{D_{(3)}}}.}
  
  \item{Y}{Response \code{\eqn{Y}}.}
  
  \item{X}{Design matrix \code{\eqn{X}}.}
  
  \item{Z}{Design matrix of Bspline approximation \code{\eqn{\lambda}}.}
  %\item{...}{ Other options for CompositeQuantile.}
}

%\author{
%Your Name, email optional.
%Maintainer: Xu Liu <liu.xu@sufe.edu.cn>
%}
\references{
  A tensor estimation approach to multivariate additive models.
}
\keyword{ High-dimensional; Sparse models; Tensor estimation; Tucker decomposition. }
\seealso{
  mam, mam_sparse_dr
}

\examples{
    D2 <- matrix(runif(50, 0.7, 1), 2, 25) # tensor with size 5*2*5
    mydata <- generateData(200, 5, 10, 5, D2)
     
    fit <- mam_sparse(mydata$Y, mydata$X)
    D3hat <- fit$Dnew
    D2hat <- TransferModalUnfoldings(D3hat,3,2,5,2,5)

}