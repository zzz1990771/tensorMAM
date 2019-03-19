\name{plotfuns}
\alias{plotfuns}
\docType{package}
\title{Plot the estimated curves from \code{tensorMam}.}
\description{
  Plot the curves fitted by \code{mam}, \code{mam_dr}, \code{mam_sparse}, and \code{mam_sparse_dr}
}

\usage{plotfuns(fit,funTrueID,true.curve=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{Object outputting from \code{mam}, \code{mam_dr}, \code{mam_sparse} or \code{mam_sparse_dr}.
}
  \item{funTrueID}{Which function to be plotted. It is a 2-vector. In MAM models, there are \eqn{s_0\times q} true functions. Thus, the first argument must be smaller than \eqn{s_0}, and the second argument must be smaller than \eqn{q}.}
  \item{true.curve}{A Logical flag. Plot both true and estimated curves if true.curve=TRUE. Plot estimated curve only if true.curve=FALSE. Default is FALSE.
  }
}

\details{
  This function gives \code{pq} functional coefficients' estimators of MAM. The singular value matrix of 
  tensor is a \eqn{r_1\times r_2\times r_3}-tensor. We choose \eqn{r_1}, \eqn{r_2}
  and \eqn{r_3} by \code{BIC} or \code{CV}.
}
\references{
  A tensor estimation approach to multivariate additive models.
}
\seealso{
  mam, mam_dr, mam_sparse, mam_sparse_dr
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
  mydata <- generateData(n, q, p, s0, D2)
  fit <- mam(mydata$Y, mydata$X)
  fit$D2 <- D2
  fit$s0 <- s0
  fit$X0 <- matrix(runif(100*p),100,p)
  plotfuns(fit, c(1,1)) 
}
