\name{mam_sparse_dr}
\alias{mam_sparse_dr-function}
\alias{mam_sparse_dr}
\docType{package}
\title{
  Fit MAM with sparsity assumption and ranks selected by \code{BIC}, \code{AIC}, \code{EBIC}, \code{CV}, or \code{GCV}.
}
\description{
  Fit a high-dimensional multivariate additive model using B-splines, with or with aparsity assumptions and ranks selected by \code{BIC}, \code{AIC}, \code{EBIC}, \code{CV}, or \code{GCV}. The multivariate sparse group \code{LASSO}, \code{MCP} or \code{SCAD}) and the coordinate descent algorithm are used to estimate
  functions for sparsity situation. The tuning parameter is selected by \code{BIC}, \code{AIC}, \code{EBIC}, \code{CV}, or \code{GCV}, which matchs the method of rank selection.
}

\usage{
mam_sparse_dr(Y, X, criteria = "BIC", ncv = 10, penalty = "LASSO", isPenColumn=TRUE, 
              r1_index = NULL, r2_index = NULL, r3_index = NULL, lambda = NULL, 
              SABC = NULL, intercept = TRUE, initMethod=NULL, nlam = 50, 
              K = 6, degr = 3, lam_min = 0.01, eps1 = 1e-4, maxstep1 = 20, 
              eps2 = 1e-4, maxstep2 = 20, gamma = 2, dfmax = NULL, alpha = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{A \eqn{n\times q} numeric matrix of responses.}
  
  \item{X}{A \eqn{n\times p} numeric design matrix for the model.}
  
  \item{criteria}{The criteria to be applied to select parameters.  Either \code{BIC} (the default), \code{AIC}, \code{EBIC}, \code{CV}, 
                or \code{GCV}.}
    
  \item{ncv}{The number of cross-validation folds.  Default is 10. If \code{criteria} is not \code{CV}, \code{ncv} is useless.}
  
  \item{penalty}{The penalty to be applied to the model. Either \code{LASSO} (the default), \code{MCP} or \code{SCAD}.}
  
  \item{isPenColumn}{A logical value indicating whether the coefficients associating with \eqn{X_j} that affects whole response 
                    \eqn{y} is penalized.  Default is \code{TRUE}. If \code{isPenColumn} is \code{TRUE}, the coefficients associating with
                    \eqn{X_j} that affects simultaneously whole response \eqn{y} is penalized for each \eqn{j\in \{1,\cdots,p\}}. 
                    If \code{isPenColumn} is \code{FALSE}, the coefficients associating with \eqn{X_j} that affects single response
                    \eqn{y_l} is penalized for each \eqn{j\in \{1,\cdots,p\}}, where \eqn{l\in \{1,\cdots,q\}}.}
                 
                 
  \item{r1_index}{A user-specified sequence of \eqn{r_1} values, where 
                  \eqn{r_1} is the first dimension of single value matrix of the  
                  tensor. Default is 
                  \code{r1_index}\eqn{=1,\cdots,\min(\lceil\log(n)\rceil,p)}.
                  }
            
  \item{r2_index}{A user-specified sequence of \eqn{r_2} values, where 
                  \eqn{r_2} is the second dimension of single value matrix of the
                  tensor. Default is \code{r2_index}\eqn{=1,\cdots},\code{max{K_index}}.
                  }
                  
  \item{r3_index}{A user-specified sequence of \eqn{r_3} values, where 
                  \eqn{r_3} is the third dimension of single value matrix of the  
                  tensor. Default is 
                  \code{r3_index}\eqn{=1,\cdots,\min(\lceil\log(n)\rceil,q)}.
                  }
                  
  \item{lambda}{A user-specified sequence of lambda values.  By default,
                a sequence of values of length \code{nlam} is computed, equally
                spaced on the log scale.}

  \item{SABC}{A user-specified list of initial coefficient matrix of \eqn{S},  
              \eqn{A}, \eqn{B}, \eqn{C}. By default,
              initial matrices are provided by random.}
              
  \item{intercept}{Should intercept(s) be fitted (default=\code{TRUE}) or set to zero (\code{FALSE})?}
  
  \item{initMethod}{One can estimate the initial tensor \eqn{D_{(3)}} as a metrix by choosing a penalty to penalize group-column wise. \code{initMethod}
                    can be \code{LASSO}, \code{MCP} or \code{SCAD}. The default is \code{NULL}}, which means no sparse initial to be provided.
                    
  \item{nlam}{The number of lambda values.  Default is 50.}
  
   \item{K}{The number of B-spline base function, that is the plus of both degrees of base  
          functioin and the number of knots. Default is \code{6}.}
		  
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
  tensor is a \eqn{r_1\times r_2\times r_3}-tensor. We choose \eqn{r_1}, \eqn{r_2} 
  and \eqn{r_3}  by \code{BIC} or \code{CV}.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use

  \item{Dnew}{Estimator of \eqn{D_{(3)}}.}
  
  \item{mu}{Estimator of intercept \eqn{\mu}.}
  
  \item{rss}{Residual sum of squares (RSS).}
  
  \item{df}{Degrees of freedom.}
  
  \item{activeX}{The active set of \eqn{X}. It is a \eqn{p} dimensional vector.}
  
  \item{activeF}{The active set of functions. It is a \eqn{q\times p} matrix.}
  
  \item{lambda}{The sequence of regularization parameter values in the path.}
  
  \item{selectedID }{The index of \code{lambda} corresponding to \code{lambda_opt}.}
  
  \item{lambda_opt }{The value of \code{lambda} with the minimum \code{BIC}, \code{AIC}, \code{EBIC}, \code{CV}, or \code{GCV} value.}
  
  \item{RSS}{The values of \code{BIC} or \code{CV}, which is a vector.}
  
  \item{rk_opt}{The optimal parametres that slected by \code{BIC}, \code{AIC}, \code{EBIC}, \code{CV}, or \code{GCV}. It is a vector with length 4, which are selected \eqn{r_1},  \eqn{r_2}, \eqn{r_3}, and \eqn{K}.}
  
  \item{Y}{Response \eqn{Y}.}
  
  \item{X}{Design matrix \eqn{X}.}
  
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

  p <- 10
  q <- 5
  s <- 5
  D2 <- matrix(runif(2*s*q, 0.7, 1), 2, s*q) # tensor with size 5*2*5
  mydata <- generateData(200, q, p, s, D2)
  
  fit <- mam_sparse_dr(mydata$Y, mydata$X)
  K <- fit$K
  D3hat <- fit$Dnew  # A q*(Kp) matrix with (p,K,q)=(5,6,5)
  D2hat <- TransferModalUnfoldings(D3hat,3,2,p,K,q)
  D1hat <- TransferModalUnfoldings(D3hat,3,1,p,K,q)
  opt <- fit$rk_opt
  which(rowSums(D1hat^2)>0)
  fit$activeX

  #Example 2

  data(breastData)
  attach(breastData)
  Y = t(dna[chrom==21,])
  Xt = t(rna[which(genechr==21),])
  n = nrow(Y)

  minX = apply(Xt,2,min)
  maxX = apply(Xt,2,max)
  X = (Xt - matrix(rep(minX,each = n),n))/matrix(rep(maxX-minX,each = n),n)
  Y = scale(Y)
  fit <- mam_sparse_dr(Y[,1:5], X[,1:10])
  K <- fit$K
  D3hat <- fit$Dnew  # A q*(Kp) matrix with (p,K,q)=(10,6,5)
  D1hat <- TransferModalUnfoldings(D3hat,3,1,10,K,5)
  opt <- fit$rk_opt
  which(rowSums(D1hat^2)>0)
  fit$activeX
  detach(breastData)
}