\name{generateData}
\alias{generateData}
\docType{package}
\title{
  Generate data from MAM model.
}
\description{
  Generate data for a high-dimensional multivariate additive model, with or without aparsity assumptions.
}

\usage{
generateData(n, q, p, s, D2, sigma2=NULL, indexF=NULL, seed_id=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Sample size.}
  
  \item{q}{The number of responses, \eqn{q\geq1}.}
  
  \item{p}{The number of covariates, \eqn{p\geq1}.}
  
  \item{s}{The true covariates associating to response, \eqn{s\geq1}.}
  
  \item{D2}{The mode of unfolding \eqn{D_{(2)}}.}
  
  \item{sigma2}{err variance. Default is 0.1.}
  
  \item{indexF}{A \eqn{q\times s} matrix. The index of significant predictors corresponding to response \eqn{y_l}. 
  Default is the matrix with each row being \eqn{(1,2,\cdots,s)}.}
  
  \item{seed_id}{A positive integer, the seed for generating the random numbers.}
}

\details{This function gives \code{pq} functional coefficients' estimators of MAM. The singular value matrix of 
  tensor is a \eqn{r_1\times r_2\times r_3}-tensor. We choose \eqn{r_1}, \eqn{r_2} 
  and \eqn{r_3}  by \code{BIC}, \code{AIC}, \code{EBIC}, \code{CV}, or \code{GCV}.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Y}{Response, a \eqn{n\times q}-matrix.}
  
  \item{X}{Design matrix, a \eqn{n\times p}-matrix.}
  
  \item{f0}{True functions.}
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
  mam_sparse
}
\examples{
  # Example 1
  
  D2 <- matrix(runif(30, 0.7, 1), 2, 15)
  mydata <- generateData(200, 3, 5, 5, D2)
    
  Y <- mydata$Y  
  X <- mydata$X
    
    
  # Example 2
  n <- 500
  p <- 10
  q <- 10
  s <- 10
  K <- 6
  s0 <- s
  r10=r20=r30=2
  S3 <- matrix(runif(r10*r20*r30,3,7),nrow = r30)
  T1 <- matrix(rnorm(s0*r10),nrow = s0)
  U1 <- qr.Q(qr(T1))
  T1 <- matrix(rnorm(K*r20),nrow = K)
  U2 <- qr.Q(qr(T1))
  T1 <- matrix(rnorm(q*r30),nrow = q)
  U3 <- qr.Q(qr(T1))
  D3 <- U3\%*\%S3\%*\%t(kronecker(U2,U1))
  D2 <- TransferModalUnfoldings(D3,3,2,s0,K,q)
  mydata <- generateData(n,q,p,s0,D2)
}
