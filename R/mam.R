
##--------------without sparsity----------------------##
mam <- function(Y,X,K=6,r1=NULL,r2=NULL,r3=NULL,SABC=NULL,intercept=TRUE,mu=NULL,degr=3,eps=1e-4,max_step=20){
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  if(is.null(r1)) r1 <- 2 
  if(is.null(r2)) r2 <- 2
  if(is.null(r3)) r3 <- 2
  if(degr>K-1) stop("K must be larger than degree+1 !")

  # initial A,B,C,S
  if(is.null(SABC)){
    set.seed(1)
    A = rbind(diag(r1), matrix(0,p-r1,r1))
    B = rbind(diag(r2), matrix(0,K-r2,r2))
    C = rbind(diag(r3), matrix(0,q-r3,r3))
    S = matrix(rnorm(r1*r2*r3),r3,r1*r2)
  }
  else{
    A = SABC$A
    B = SABC$B
    C = SABC$C
    S = SABC$S
  }
  if(intercept & is.null(mu)) mu = as.vector(rep(0,q))
  
  Z = bsbasefun(X,K,degr)
  fit = Estimation(Y,Z,A,B,C,S,intercept,mu,eps,max_step)
  return(list(Dnew=fit$Dnew, 
              rss=fit$likhd,
              mu = fit$mu,
              Y = Y,
              X = X,
              Z = Z,
              degr = degr,
              K = K
              )
         )
}