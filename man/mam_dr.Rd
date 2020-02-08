\name{mam_dr}
\alias{mam_dr-function}
\alias{mam_dr}
\docType{package}
\title{
  Fit MAM without sparsity assumption, and with ranks selected by \code{BIC}, \code{AIC}, \code{EBIC}, \code{CV}, or \code{GCV}.
}
\description{
  Fit a low-dimensional multivariate additive model using B-splines, without aparsity assumptions, and with ranks \eqn{r_1, r_2, r_3} selected by \code{BIC}, \code{AIC}, \code{EBIC}, \code{CV}, or \code{GCV}.
} 

\usage{mam_dr(Y, X, criteria = "BIC", ncv = 10, r1_index = NULL, r2_index = NULL, 
       r3_index = NULL, SABC = NULL, intercept = TRUE, K = 6, degr = 3, 
       eps = 1e-4, max_step = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{A \eqn{n\times q} numeric matrix of responses.}
  
  \item{X}{A \eqn{n\times p} numeric design matrix for the model.}

  \item{criteria}{The criteria to be applied to select parameters.  Either \code{BIC} (the default), \code{AIC}, \code{EBIC}, \code{CV}, 
                or \code{GCV}.}
    
  \item{ncv}{The number of cross-validation folds.  Default is 10. If \code{criteria} is not  
             \code{"CV"}, \code{ncv} is useless.}
 
    
  \item{r1_index}{A user-specified sequence of \eqn{r_1} values, where 
                   \eqn{r_1} is the first dimension of single value matrix of the 
                   tensor. Default is 
                   \code{r1_index}\eqn{=1,\cdots,\min(\log(n)\rceil,p)}.
                   }
   
  \item{r2_index}{A user-specified sequence of \eqn{r_2} values, where 
                  \eqn{r_2} is the second dimension of single value matrix of the 
                  tensor. Default is \code{r2_index}\eqn{=1,\cdots},\code{max{K_index}}.
                  }
   
  \item{r3_index}{A user-specified sequence of \eqn{r_3} values, where 
                  \eqn{r_3} is the third dimension of single value matrix of the  
                  tensor. Default is 
                  \code{r3_index}\eqn{=1,\cdots,\min(\log(n)\rceil,q)}.
                  }
                   
  \item{SABC}{A user-specified list of initial coefficient matrix of \eqn{S},
              \eqn{A}, \eqn{B}, \eqn{C}, which is a list with values \eqn{S},
              \eqn{A}, \eqn{B}, \eqn{C}. By default,
              initial matrices are provided by random.}
        
  \item{intercept}{Should intercept(s) be fitted (default=\code{TRUE}) or set to zero (\code{FALSE})?}
  
  \item{K}{The number of B-spline base function, that is the plus of both degrees of base  
          functioin and the number of knots. Default is \code{6}.}
  
  \item{degr}{ the number of knots of B-spline base function. Default is \code{3}.}
  
  \item{eps}{Convergence threshhold.  The algorithm iterates until the
             relative change in any coefficient is less than \code{eps}.  Default
             is \code{1e-4}.}
             
  \item{max_step}{Maximum number of iterations.  Default is 20.}
}

\details{This function gives \code{pq} functional coefficients' estimators of MAM. The singular value matrix of 
  tensor is a \eqn{r_1\times r_2\times r_3}-tensor. We choose \eqn{r_1}, \eqn{r_2} 
  and \eqn{r_3}  by \code{BIC} or \code{CV}.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Dnew}{Estimator of \eqn{D_{(3)}}.}
  
  \item{mu}{Estimator of intercept \eqn{\mu}.}
  
  \item{rss}{Residual sum of squares (RSS).}
  
  \item{rk_opt}{The optimal parametres that slected by \code{BIC}, \code{AIC}, \code{EBIC}, \code{CV}, or \code{GCV}. It is a vector with length 4, which are selected \eqn{r_1},  \eqn{r_2}, \eqn{r_3}, and \eqn{K}.}
  
  \item{selected}{Which \eqn{\lambda} is selection.}
  
  \item{Y}{Response \eqn{Y}.}
  
  \item{X}{Design matrix \eqn{X}.}
  
  \item{Z}{Design matrix of Bspline approximation.}
}

%\author{
%Your Name, email optional.
%Maintainer: Xu Liu <liu.xu@sufe.edu.cn>
%}
\references{A tensor estimation approach to multivariate additive models.
}
\seealso{
  mam, mam_sparse_dr
}
\examples{ 
  p <- 5
  q <- 5
  D2 <- matrix(runif(2*p*q, 0.7, 1), 2, p*q) # tensor with size 5*2*5
  mydata <- generateData(200, q, p, p, D2)
  
  fit <- mam_dr(mydata$Y, mydata$X)
  K <- fit$K
  D3hat <- fit$Dnew  # A q*(Kp) matrix with (p,K,q)=(5,6,5)
  D2hat <- TransferModalUnfoldings(D3hat,3,2,p,K,q)
  opt <- fit$rk_opt
}