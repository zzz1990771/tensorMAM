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
  \item{fit}{??????\code{mam} or \code{mam_sparse} object.
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
\dontrun{
    D2 <- matrix(runif(50, 0.7, 1), 2, 25)
    mydata <- generateData(200, 5, 5, 5, D2)
    fit <- mam(mydata$Y, mydata$X)
    plotfuns(fit, c(1,1)) 
}
}
