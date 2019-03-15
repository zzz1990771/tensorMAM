# list of functions which are used by R codes.
assig <- function(n_args){
  cargs <- vector("list", length(n_args))
  for(i in 1:length(n_args)) cargs[[i]] <- 1:n_args[i]
  t(expand.grid(cargs))
}

##--------------produce the B-spline functions----------------------##
bsbasefun1 <- function(X,K,degr){
  n = dim(X)[1]
  p = dim(X)[2]
  nk = K - degr - 1
  u.k = seq(0, 1, length=nk+2)[-c(1,nk+2)]
  BS = NULL
  for(j in 1:p){
    Knots = as.numeric(quantile(X[,j], u.k))  
    BS0 = bs(X[,j], knots=Knots, intercept=TRUE, degree=degr)
    BS = cbind(BS,BS0)
  }
  id = seq(1,p*K,K)
  Z = NULL
  for(j in 1:K){
    Z = cbind(Z,BS[,id+j-1])
  }
  return(Z)
}

bsbasefun <- function(X,K,degr){
  n = dim(X)[1]
  p = dim(X)[2]
  nk = K - degr
  u.k = seq(0, 1, length=nk+2)[-c(1,nk+2)]
  BS = NULL
  for(j in 1:p){
    Knots = as.numeric(quantile(X[,j], u.k))  
    BS0 = bs(X[,j], knots=Knots, intercept=TRUE, degree=degr)
    BS = cbind(BS,BS0[,-1])
  }
  BS = scale(BS,center = T, scale = F)
  id = seq(1,p*K,K)
  Z = NULL
  for(j in 1:K){
    Z = cbind(Z,BS[,id+j-1])
  }
  return(Z)
}

##--------------produce p*q estimated functions----------------------##
trans <- function(X,D3,p,q,K,degr,s0){
  Z = bsbasefun(X,K,degr)
  id = seq(1,p*K,p)
  fjl = assig(c(s0,q))
  funhat = NULL
  for(j0 in 1:(q*s0)){
    qj1 = fjl[1,j0]
    qj = fjl[2,j0]
    funhat = cbind(funhat,Z[,id+qj1-1]%*%D3[qj,id+qj1-1])
  }
  return(funhat)
}

truefuns <- function(n,q,p,s,D2,seed_id){
  set.seed(seed_id)
  X <- matrix(runif(n*p),n,p)
  X1 <- X[,1:s]
  basefuns1 <- sin(2*pi*X1)
  basefuns2 <- cos(pi*X1)
  f0 <- matrix(rep(basefuns1,q),nrow=n)*matrix(rep(D2[1,],each=n),n) + matrix(rep(basefuns2,q),nrow=n)*matrix(rep(D2[2,],each=n),n)
  return(list(X=X,f0=f0))
  
}
##--------------generate data----------------------##
generateData1 <- function(n,q,p,s,D2,seed_id,rho=0.0,S = NULL){
  set.seed(seed_id)
  if(is.null(S)) {
    S = diag(p)
    for(j in 1:p){
      for(i in j:p) S[i,j]= rho^(i-j)
    }
    S = S + t(S) - diag(p)
  }
  Sroot <- chol(S)
  X <- matrix(runif(n*p), nrow = n)%*%Sroot
  X1 <- X[,1:s]
  basefuns1 <- sin(2*pi*X1)
  basefuns2 <- cos(pi*X1)
  
  fl <- basefuns1%*%matrix(D2[1,],nrow=s)+basefuns2%*%matrix(D2[2,],nrow=s)
  f0 <- matrix(rep(basefuns1,q),nrow=n)*matrix(rep(D2[1,],each=n),n) + matrix(rep(basefuns2,q),nrow=n)*matrix(rep(D2[2,],each=n),n)
  
  eps <- matrix(rnorm(n*q),n,q)
  Y <- fl + eps*0.1
  return(list(Y=Y,X=X,f0=f0))
}

##--------------generate data----------------------##
generateData <- function(n,q,p,s,D2, sigma2=NULL, t=NULL,seed_id=1e4){
  if(is.null(t)) t = 0
  if(is.null(sigma2)) sigma2 = 0.1
  set.seed(seed_id)
  rho = 0.0
  W <- matrix(runif(n*p), nrow = n)
  U <- runif(n)
  X <- (W+t*matrix(rep(U,p),nrow = n))/(1+t)
  X1 <- X[,1:s]
  basefuns1 <- sin(2*pi*X1)
  basefuns2 <- cos(pi*X1)
  
  fl <- basefuns1%*%matrix(D2[1,],nrow=s)+basefuns2%*%matrix(D2[2,],nrow=s)
  f0 <- matrix(rep(basefuns1,q),nrow=n)*matrix(rep(D2[1,],each=n),n) + matrix(rep(basefuns2,q),nrow=n)*matrix(rep(D2[2,],each=n),n)
  
  Sigma = matrix(rho,q,q)+diag(1-rho,q,q)
  A = chol(Sigma) #t(A)*A=Sigma
  eps <- matrix(rnorm(n*q),n,q)
  Y <- fl + sigma2*eps%*%A
  #Y <- fl + eps
  return(list(Y=Y,X=X,f0=f0))
}

##--------------generate data----------------------##
generateData2 <- function(n,q,p,s,D3,K,degr,seed_id, t=NULL, sigma2=NULL){
  if(is.null(t)) t = 0
  if(is.null(sigma2)) sigma2 = 0.1
  set.seed(seed_id)
  W <- matrix(runif(n*p), nrow = n)
  U <- runif(n)
  X <- (W+t*matrix(rep(U,p),nrow = n))/(1+t)
  X1 <- X[,1:s]
  
  Z = bsbasefun(X1,K,degr)
  fl = Z%*%t(D3)
  f0 = fl
  eps <- matrix(rnorm(n*q),n,q)
  Y <- fl + sigma2 * eps
  return(list(Y=Y,X=X,f0=f0))
}

##--------------plot curve of function f_{jl} ----------------------##
plotfuns <- function(fit,funTrueID){
  # funTrueID = c(j,l) is the index of the f_{jl}th function
  n = fit$n
  p = fit$p
  q = fit$q
  s0 = fit$s0
  D2 = fit$D2
  X = fit$X0
  funhat = fit$funhat
  n0 = dim(X)[1]
  
  qj1 = funTrueID[1]
  qj = funTrueID[2]
  j0 = s0*(qj-1) + qj1
  if(qj1>s0) stop("j must not be larger than s !")
  if(qj>q) stop("l must not be larger than q !")
  
  X1 <- X[,1:s0]
  basefuns1 <- sin(2*pi*X1)
  basefuns2 <- cos(pi*X1)
  f0 <- matrix(rep(basefuns1,q),nrow=n0)*matrix(rep(D2[1,],each=n0),n0)+matrix(rep(basefuns2,q),nrow=n0)*matrix(rep(D2[2,],each=n0),n0)
  xs = sort(X[,qj1],index.return = T)
  w = xs$x
  iw = xs$ix
  if(qj1>s0)  f11 = rep(0,n)
  if(qj1<=s0) f11 = f0[,j0]
  
  f11hat = funhat[,j0]
  plot(w,f11[iw],type = "l",col = "blue",ylim=c(min(f11[iw])-0.2,max(f11[iw])+0.2))
  lines(w,f11hat[iw],col="red")
}


Imse <- function(fit){
  n = fit$n
  p = fit$p
  q = fit$q
  s0 = fit$s0
  D2 = fit$D2
  X = fit$X0
  funhat = fit$funhat
  
  n0 = dim(X)[1]
  X1 <- X[,1:s0]
  basefuns1 <- sin(2*pi*X1)
  basefuns2 <- cos(pi*X1)
  f0 <- matrix(rep(basefuns1,q),nrow=n0)*matrix(rep(D2[1,],each=n0),n0)+matrix(rep(basefuns2,q),nrow=n0)*matrix(rep(D2[2,],each=n0),n0)
  return(colSums((f0-funhat)^2)/n0)
}

plotfuns1 <- function(fit,funTrueID){
  # funTrueID = c(j,l) is the index of the f_{jl}th function
  n = fit$n
  p = fit$p
  q = fit$q
  s0 = fit$s0
  D2 = fit$D2
  X = fit$X0
  funhat = fit$funhat
  
  n0 = dim(X)[1]
  
  # Let true functions is a matrix F with dimension p*q
  # Then qj is the column number of F
  # qj1 is the row number of F
  # that is, (qj1,qj) element of F
  j0 = funTrueID
  qj = ceiling(j0/s0)
  qj1 = j0%%s0
  
  if(qj1==0) qj1=s0
  if(qj1>s0) stop("j must not be larger than s !")
  if(qj>q) stop("l must not be larger than q !")
  cat("(j,l) = ",qj1,qj,"\n")
  
  
  X1 <- X[,1:s0]
  basefuns1 <- sin(2*pi*X1)
  basefuns2 <- cos(pi*X1)
  f0 <- matrix(rep(basefuns1,q),nrow=n0)*matrix(rep(D2[1,],each=n0),n0)+matrix(rep(basefuns2,q),nrow=n0)*matrix(rep(D2[2,],each=n0),n0)
  xs = sort(X[,qj1],index.return = T)
  w = xs$x
  iw = xs$ix
  if(qj1>s0)  f11 = rep(0,n)
  if(qj1<=s0) f11 = f0[,j0]
  
  f11hat = funhat[,j0]
  plot(w,f11[iw],type = "l",col = "blue",ylim=c(min(f11[iw])-0.1,max(f11[iw])+0.1))
  lines(w,f11hat[iw],col="red")
}

plotfuns2 <- function(fit,funTrueID){
  # funTrueID = c(j,l) is the index of the f_{jl}th function
  n = fit$n
  p = fit$p
  q = fit$q
  s0 = fit$s0
  D2 = fit$D2
  X = fit$X0
  K = fit$K
  funhat = fit$funhat
  
  degr = 3
  n0 = dim(X)[1]
  
  # Let true functions is a matrix F with dimension p*q
  # Then qj is the column number of F
  # qj1 is the row number of F
  # that is, (qj1,qj) element of F
  j0 = funTrueID
  qj = ceiling(j0/s0)
  qj1 = j0%%s0
  
  if(qj1==0) qj1=s0
  if(qj1>s0) stop("j must not be larger than s !")
  if(qj>q) stop("l must not be larger than q !")
  cat("(j,l) = ",qj1,qj,"\n")
  
  
  X1 <- X[,1:s0]
  D3 = TransferModalUnfoldings(D2,2,3,s0,K,q)
  f0 = trans(X,D3,p,q,K,degr,s0)
  xs = sort(X[,qj1],index.return = T)
  w = xs$x
  iw = xs$ix
  if(qj1>s0)  f11 = rep(0,n)
  if(qj1<=s0) f11 = f0[,j0]
  
  f11hat = funhat[,j0]
  plot(w,f11[iw],type = "l",col = "blue",ylim=c(min(f11[iw])-0.1,max(f11[iw])+0.1))
  lines(w,f11hat[iw],col="red")
}

Imse2 <- function(fit){
  n = fit$n
  p = fit$p
  q = fit$q
  s0 = fit$s0
  D2 = fit$D2
  X = fit$X0
  K = fit$K
  degr = fit$degr
  funhat = fit$funhat
  
  n0 = dim(X)[1]
  X1 <- X[,1:s0]
  D3 = TransferModalUnfoldings(D2,2,3,s0,K,q)
  f0 = trans(X1,D3,s0,q,K,degr,s0)
  return(colSums((f0-funhat)^2)/n0)
}


