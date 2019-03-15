\name{plotfuns}
\alias{plotfuns}
\docType{package}
\title{Plot the estimated curves from \codes{tensorMam}.}
\description{
  Plot the curves fit by \code{mam}, \code{mam_dr}, \code{mam_sparse}, and \code{mam_sparse_dr}
}

\usage{plotfuns(fit,funTrueID)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{\code{mam} or \code{mam_sparse} object.
}
  \item{funTrueID}{Which function be plotted.}
}

\details{
  This function gives \code{pq} functional coefficients' estimators of MAM. The singular value matrix of 
  tensor is a \code{\eqn{r_1\times r_2\times r_3}}-tensor. We choose \code{\eqn{r_1}}, \code{\eqn{r_2}} 
  and \code{\eqn{r_3}}  by \code{BIC} or \code{CV}.
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
  mam
}
\examples{
  n <- 200
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
  mydata <- generateData(n, q, p, K, D2)
  fit <- mam(mydata$Y, mydata$X)
  fit$D2 <- D2
  fit$s0 <- s0
  fit$X0 <- matrix(runif(100*p),100,p)
  plotfuns(fit, c(1,1)) 
}
