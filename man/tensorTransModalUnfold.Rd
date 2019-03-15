\name{TransferModalUnfoldings}
\alias{TransferModalUnfoldings-function}
\alias{TransferModalUnfoldings}
\docType{package}
\title{Transfer a tensor's modal unfoldings to another.}
\description{
  Transfer a tensor's modal unfoldings to another.
}

\usage{
tensorTransModeUnfold(T, d1, d2 , r1, r2, r3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{
  A mode-d1-unfolding of a tensor with size \code{\eqn{r_1\times r_2\times r_3}}, input unfolding
}
  \item{d1}{
  An integer, the mode of unfolding \code{\eqn{T_{(d_1)}}}
}
  \item{d2}{
   An integer, the mode of output unfolding \code{\eqn{T_{(d_2)}}}
}
  \item{r1}{
  The fist dimension of tensor 
}
  \item{r2}{
  The second dimension of tensor
}
  \item{r3}{
  The third dimension of tensor 
}
}

\details{
  This function transfers an input mode-d1-unfolding \code{\eqn{T_{(d_1)}}} 
  to mode-d2-unfolding \code{\eqn{T_{(d_2)}}}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{D}{the output mode-d2-unfolding, \code{\eqn{T_{(d_2)}}}}
}

%\author{
%Your Name, email optional.
%Maintainer: Xu Liu <liu.xu@sufe.edu.cn>
%}
\references{
  A tensor estimation approach to multivariate additive models.
}
\keyword{ High-dimensional, Sparse models; Tensor estimation; Tucker decomposition. }
\examples{
  T <- matrix(1:24,nrow = 4)
  TransferModalUnfoldings(T,1,2,4,3,2)
}
