\name{mam}
\alias{mam-function}
\alias{mam}
\docType{package}
\title{
  Fit MAM without sparsity assumption and with fixed ranks.
}
\description{
  Fit a low-dimensional multivariate additive model using B-splines, without aparsity assumptions, and given ranks given ranks \eqn{r_1, r_2, r_3}.
}

\usage{mam(Y, X, K = 6, r1 = NULL, r2 = NULL, r3 = NULL, SABC = NULL, degr = 3,
       eps = 1e-4, max_step = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{A \eqn{n\times q} numeric matrix of responses.}
  
  \item{X}{A \eqn{n\times p} numeric design matrix for the model.}
  
  \item{K}{The number of B-spline base function, that is the plus of both degrees of base functioin and the number of knots. Default is \code{6}.}
  
  \item{degr}{ the number of knots of B-spline base function. Default is \code{3}.}
  
  \item{r1}{The first dimension of single value matrix of the tensor. Default is 2.}
  
  \item{r2}{The second dimension of single value matrix of the tensor. Default is 2.}
  
  \item{r3}{The third dimension of single value matrix of the tensor. Default is 2.}
  
  \item{SABC}{A user-specified list of initial coefficient matrix of \eqn{S},
              \eqn{A}, \eqn{B}, \eqn{C}. By default,
              initial matrices are provided by random.}

  \item{eps}{Convergence threshhold.  The algorithm iterates until the
    relative change in any coefficient is less than \code{eps}.  Default
    is \code{1e-4}.}
    
  \item{max_step}{Maximum number of iterations.  Default is 20.}
}

\details{This function gives \code{pq} functional coefficients' estimators of MAM. The singular value matrix of 
  tensor is a \code{\eqn{r_1\times r_2\times r_3}}-tensor. We choose \eqn{r_1}, \eqn{r_2} 
  and \eqn{r_3}  by \code{BIC} or \code{CV}.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Dnew}{Estimator of \eqn{D_{(3)}}.}
  
  \item{rss }{Residual sum of squares (RSS).}
  
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
\keyword{High-dimensional; Sparse models; Tensor estimation; Tucker decomposition. }
\seealso{
  mam_sparse
}
\examples{ 
  D2 <- matrix(runif(50, 0.7, 1), 2, 25) # tensor with size 5*2*5
  mydata <- generateData(200, 5, 5, 5, D2)
    
  fit <- mam(mydata$Y, mydata$X)
  D3hat <- fit$Dnew
  D2hat <- TransferModalUnfoldings(D3hat,3,2,5,2,5)
}