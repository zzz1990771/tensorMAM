
##--------------main by BIC without sparsity----------------------##
mam_dr <- function(Y,X,criteria="BIC",ncv=10,r1_index=NULL,r2_index=NULL,r3_index=NULL,SABC=NULL,
                   intercept=TRUE,K=6,degr=3,eps=1e-4,max_step=20){

  if(degr>K-1) stop("K must be larger than degree+1 !")
  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)
  if(is.null(r1_index)) r1_index = 1:min(ceiling(log(n)),p)
  if(is.null(r2_index)) r2_index = 1:min(K)
  if(is.null(r3_index)) r3_index = 1:min(ceiling(log(n)),q)
  
  #---------------- The selection by BIC  ---------------------#  
  if(is.null(SABC)){
    set.seed(1)
    r1_max = max(r1_index) 
    r2_max = max(r2_index) 
    r3_max = max(r3_index) 
    A = rbind(diag(r1_max), matrix(0,p-r1_max,r1_max))
    B = rbind(diag(r2_max), matrix(0,K-r2_max,r2_max))
    C = rbind(diag(r3_max), matrix(0,q-r3_max,r3_max))
    S = matrix(rnorm(r1_max*r2_max*r3_max),r3_max,r1_max*r2_max)
  }
  else{
    A = SABC$A
    B = SABC$B
    C = SABC$C
    S = SABC$S
  }
  if((max(r1_index)>dim(A)[2])|(max(r2_index)>dim(B)[2])|(max(r3_index)>dim(C)[2]))
    stop("maximum number of index sequence of r1, r2, and r3 must not be larger than A, B, and C, respectively !")

  if(criteria=="CV") fit_dr = mam_cv(Y,X,ncv,r1_index,r2_index,r3_index,S,A,B,C,intercept,K,degr,eps,max_step)
  else fit_dr = mam_bic(Y,X,criteria,r1_index,r2_index,r3_index,S,A,B,C,intercept,K,degr,eps,max_step)
  
  return(fit_dr)
}