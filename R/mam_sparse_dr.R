
##--------------Estimation with Penalty by CV----------------------##
mam_sparse_dr <- 
  function(Y,X,method="BIC",ncv=10,penalty="LASSO",isPenColumn=TRUE,K_index=NULL,r1_index=NULL,r2_index=NULL,r3_index=NULL,lambda=NULL,
           SABC=NULL,intercept=TRUE,mu=NULL,nlam=50,degr=3,lam_min=0.01,eps1=1e-4,maxstep1=20,eps2=1e-4,maxstep2=20,gamma=2,dfmax=NULL,alpha=1){
    n <- dim(Y)[1]
    q <- dim(Y)[2]
    p <- dim(X)[2]
    K1 <- 7
    if(degr>min(6,K1-1)-1) stop("K must be larger than degree+1 !")
    if(is.null(K_index)) K_index = min(6,K1-1):max(8,K1+1)
    if(is.null(r1_index)) r1_index = 1:min(floor(log(n)),p)
    if(is.null(r2_index)) r2_index = 1:min(K_index)
    if(is.null(r3_index)) r3_index = 1:min(floor(log(n)),q)
    if (penalty == "LASSO") pen <- 1
    if (penalty == "MCP")   pen <- 2 
    if (penalty=="SCAD"){    
      gamma <- 3
      pen <- 3
    }  
    if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty")
    if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty")
    
    if (is.null(dfmax)) dfmax = p + 1
    # initial A,B,C,S
    if(is.null(SABC)){
      set.seed(1)
      r1_max = max(r1_index) 
      r2_max = max(r2_index) 
      r3_max = max(r3_index) 
      K_max = max(K_index)
      A = rbind(diag(r1_max), matrix(0,p-r1_max,r1_max))
      B = rbind(diag(r2_max), matrix(0,K_max-r2_max,r2_max))
      C = rbind(diag(r3_max), matrix(0,q-r3_max,r3_max))
      S = matrix(rnorm(r1_max*r2_max*r3_max),r3_max,r1_max*r2_max)
    }
    else{
      A = SABC$A
      B = SABC$B
      C = SABC$C
      S = SABC$S
    }
    if(is.null(mu)) mu = as.vector(rep(0,q))
    if(!intercept) mu = as.vector(rep(0,q))
    
    if (is.null(lambda)) {
      is_setlam = 1
      if (nlam < 1||is.null(nlam)) stop("nlambda must be at least 1")
      if (n<=p) lam_min = 1e-1
      setlam = c(1,lam_min,alpha,nlam)
      Z = bsbasefun(X,max(K_index),degr)
      lambda = setuplambda(Y,Z,A,B,C,S,nlam,setlam)
    }
    else {
      is_setlam = 0
      nlam = length(lambda)
      setlam = c(1,lam_min,alpha,nlam)
    }
    #---------------- The selection by CV  ---------------------#  
    if((max(r1_index)>dim(A)[2])|(max(r2_index)>dim(B)[2])|(max(r3_index)>dim(C)[2]))
      stop("maximum number of index sequence of r1, r2, and r3 must not be larger than A, B, and C, respectively !")
    
    if(method=="CV") fit_dr = mam_sparse_cv(Y,X,ncv,K_index,r1_index,r2_index,r3_index,pen,isPenColumn,lambda,A,B,C,S,
                                            intercept,mu,nlam,degr,lam_min,eps1,maxstep1,eps2,maxstep2,gamma,dfmax,alpha)
    else fit_dr = mam_sparse_bic(Y,X,method,K_index,r1_index,r2_index,r3_index,pen,isPenColumn,lambda,A,B,C,S,
                                 intercept,mu,nlam,degr,lam_min,eps1,maxstep1,eps2,maxstep2,gamma,dfmax,alpha)

    
    return(fit_dr)
  }