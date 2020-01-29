\name{mam_sparse}
\alias{mam_sparse-function}
\alias{mam_sparse}
\docType{package}
\title{
  Fit MAM with sparsity assumption and fixed ranks.
}
\description{
  Fit a high-dimensional multivariate additive model using B-splines, with or without aparsity assumptions, and given ranks \eqn{r_1, r_2, r_3}. The multivariate sparse group \code{LASSO}, \code{MCP} or \code{SCAD}) and the coordinate descent algorithm are used to estimate functions for sparsity situation.
}

\usage{
mam_sparse(Y, X, method="BIC",ncv=10, K = 6, r1 = NULL, r2 = NULL, r3 = NULL,
           penalty="LASSO", isPenColumn=TRUE, lambda = NULL, SABC = NULL, intercept = TRUE,
           initMethod="LASSO", degr = 3, nlam = 20, lam_min = 1e-3, eps1 = 1e-4, maxstep1 = 20, 
           eps2 = 1e-4, maxstep2 = 20, gamma = 2, dfmax = NULL, alpha = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{A \eqn{n\times q} numeric matrix of responses.}
  
  \item{X}{A \eqn{n\times q} numeric design matrix for the model.}
  
  \item{method}{The method to be applied to select parameters.  Either \code{BIC}
                (the default), \code{AIC}, \code{EBIC}, \code{CV}, or \code{GCV}.}
    
  \item{ncv}{The number of cross-validation folds.  Default is 10. If \code{method} is not \code{CV}, \code{ncv} is useless.}  
  
  \item{penalty}{The penalty to be applied to the model. Either \code{LASSO} (the default),  
                \code{MCP} or \code{SCAD}.}
 
  \item{isPenColumn}{A logical value indicating whether the coefficients associating with \eqn{X_j} that affects whole response 
                    \eqn{y} is penalized.  Default is \code{TRUE}. If \code{isPenColumn} is \code{TRUE}, the coefficients associating with
                    \eqn{X_j} that affects simultaneously whole response \eqn{y} is penalized for each \eqn{j\in \{1,\cdots,p\}}. 
                    If \code{isPenColumn} is \code{FALSE}, the coefficients associating with \eqn{X_j} that affects single response
                    \eqn{y_l} is penalized for each \eqn{j\in \{1,\cdots,p\}}, where \eqn{l\in \{1,\cdots,q\}}.} 
 
  \item{K}{The number of B-spline base function, that is the plus of both degrees of base  
          functioin and the number of knots. Default is \code{6}.}
          
  \item{degr}{The number of knots of B-spline base function. Default is \code{3}.}
  
  \item{r1}{The first dimension of single value matrix of the tensor. Default is 2.}
  
  \item{r2}{The second dimension of single value matrix of the tensor. Default is 2.}
  
  \item{r3}{The third dimension of single value matrix of the tensor. Default is 2.}
  
  \item{lambda}{A user-specified sequence of lambda values.  By default,
        a sequence of values of length \code{nlam} is computed, equally
        spaced on the log scale.}
        
  \item{SABC}{A user-specified list of initial coefficient matrix of \eqn{S},
              \eqn{A}, \eqn{B}, \eqn{C}. By default,
              initial matrices are provided by random.}
   
  \item{intercept}{Should intercept(s) be fitted (default=\code{TRUE}) or set to zero (\code{FALSE})?}
  
  \item{initMethod}{One can estimate the initial tensor \eqn{D_{(3)}} as a metrix by choosing a penalty to penalize group-column wise. \code{initMethod}
                    can be \code{LASSO}, \code{MCP} or \code{SCAD}. The default is \code{LASSO}}
  
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
  tensor is a \eqn{r_1\times r_2\times r_3}-tensor. \eqn{r_1}, \eqn{r_2} 
  and \eqn{r_3}  are fixed.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{betapath}{Solution path of \eqn{\beta}.}
  
  \item{rss}{Residual sum of squares (RSS).}
  
  \item{df}{Degrees of freedom.}
  
  \item{lambda}{The sequence of regularization parameter values in the path.}
  
  \item{lambda_opt}{The value of \code{lambda} with the minimum
                    \code{BIC} value.}
  
  \item{selectedID}{The index of \code{lambda} corresponding to
                    \code{lambda_opt}.}
  
  \item{activeX}{The active set of \eqn{X}. It is a \eqn{p} dimensional vector.}
  
  \item{activeF}{The active set of functions. It is a \eqn{q\times p} matrix.}
  
  \item{Dnew}{Estimator of \eqn{D_{(3)}}.}
  
  \item{mu}{Estimator of intercept \eqn{\mu}.}
  
  \item{Y}{Response \eqn{Y}.}
  
  \item{X}{Design matrix \eqn{X}.}
  
  \item{Z}{Design matrix of Bspline approximation \eqn{\lambda}.}
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
  p <- 10
  q <- 5
  s <- 5
  D2 <- matrix(runif(2*s*q, 0.7, 1), 2, s*q) # tensor with size 5*2*5
  mydata <- generateData(200, q, p, s, D2)
   
  fit <- mam_sparse(mydata$Y, mydata$X)
  K <- fit$K
  D3hat <- fit$Dnew  # A q*(Kp) matrix with (p,K,q)=(10,6,5)
  D2hat <- TransferModalUnfoldings(D3hat,3,2,p,K,q)
  D1hat <- TransferModalUnfoldings(D3hat,3,1,p,K,q)
  which(rowSums(D1hat^2)>0)
  fit$activeX
}