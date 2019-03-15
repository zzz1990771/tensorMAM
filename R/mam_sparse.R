mam_sparse <- 
  function(Y,X,K=6,r1=NULL,r2=NULL,r3=NULL,penalty="LASSO",lambda=NULL,SABC=NULL,degr=3,nlam=20, 
           lam_min=1e-3, eps1=1e-4,maxstep1=20,eps2=1e-4,maxstep2=20,gamma=2,dfmax=NULL,alpha=1){
    n <- dim(Y)[1]
    q <- dim(Y)[2]
    p <- dim(X)[2]
    if(degr>K-1) stop("K must be larger than degree+1 !")
    if(is.null(r1)) r1 <- 2 
    if(is.null(r2)) r2 <- 2
    if(is.null(r3)) r3 <- 2
    if (penalty == "LASSO") pen <- 1
    if (penalty == "MCP")   pen <- 2 
    if (penalty=="SCAD"){    
      gamma <- 3
      pen <- 3;
    }  
    if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty")
    if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty")
    
    if (is.null(dfmax)) dfmax = p + 1
    if (is.null(lambda)) {
      is_setlam = 1
      if (nlam < 1) stop("nlambda must be at least 1")
      if (n<=p) lam_min = 1e-2
      lambda = 0
    }
    else {
      is_setlam = 0
      nlam = length(lambda)
    }
    setlam = c(1,lam_min,alpha,nlam)
    
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
    Z = bsbasefun(X,K,degr)
    fit = Estimation_penalty(Y,Z,A,B,C,S,lambda,alpha, gamma, pen, dfmax, 
                             eps1, eps2, maxstep1, maxstep2,is_setlam,setlam) 
    
    if(nlam>1){
      df = fit$df*r1
      bic = 2*log(fit$likhd) + log(n)*df/n
      selected = which.min(bic)
      lambda_opt = fit$lambda[selected]
      activeA = fit$betapath[,selected]
      fit_opt = Estimation_penalty(Y,Z,as.matrix(A),as.matrix(B),as.matrix(C),as.matrix(S),lambda[1:selected],alpha, gamma, pen, dfmax,
                                   eps1, eps2, maxstep1, maxstep2,is_setlam,setlam)
      
      Dnew = fit_opt$Dnew
    }
    else{ 
      selected = 1
      Dnew = fit$Dnew
      lambda_opt = lambda
      activeA = fit$betapath
    }
    return(list(betapath=fit$betapath, 
                rss=fit$likhd[selected],
                df = fit$df,
                lambda = fit$lambda,
                lambda_opt=lambda_opt,
                selectedID = selected,
                activeA = activeA,
                Dnew=Dnew,
                Y = Y,
                X = X,
                Z = Z
                )
           )
  }