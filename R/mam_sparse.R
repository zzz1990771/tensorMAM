mam_sparse <- 
  function(Y,X,K=6,r1=NULL,r2=NULL,r3=NULL,method="BIC",ncv=10,penalty="LASSO",isPenColumn=TRUE,lambda=NULL,SABC=NULL,
           intercept=TRUE,mu=NULL,degr=3,nlam=20,lam_min=1e-3, eps1=1e-4,maxstep1=20,eps2=1e-4,maxstep2=20,gamma=2,dfmax=NULL,alpha=1){
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
    if(is.null(mu)) mu = as.vector(rep(0,q))
    if (is.null(lambda)) {
      is_setlam = 1
      if (nlam < 1||is.null(nlam)) stop("nlambda must be at least 1")
      if (n<=p) lam_min = 1e-1
      setlam = c(1,lam_min,alpha,nlam)
      Z = bsbasefun(X,K,degr)
      lambda = setuplambda(Y,Z,A,B,C,S,nlam,setlam)
    }
    else {
      is_setlam = 0
      nlam = length(lambda)
      setlam = c(1,lam_min,alpha,nlam)
    }
    #---------------- The selection by BIC or CV  ---------------------# 
    Z = bsbasefun(X,K,degr)
    if(method=="BIC"){
      if(isPenColumn){
        fit = EstPenColumn(Y,Z,as.matrix(A),as.matrix(B),as.matrix(C),as.matrix(S),
                           intercept,mu,lambda, alpha, gamma, pen, dfmax, eps1, eps2, maxstep1, maxstep2)
        df = fit$df*r1
        loglikelih = (n*q)*log(fit$likhd/(n*q))
        bic <- switch (method,
                       BIC = loglikelih + log(n*q)*df,
                       AIC = loglikelih + 2*df,
                       GCV = fit$likhd*(n*q)/(n*q-df)^2,
                       EBIC = loglikelih + log(n*q)*df + 2*(lgamma(q*p*(p+1)/2+1) 
                                                            - lgamma(df+1) - lgamma(q*p*(p+1)/2-df+1))
        )
      }
      else{
        fit = EstPenSingle(Y,Z,as.matrix(A),as.matrix(B),as.matrix(C),as.matrix(S),
                           intercept,mu,lambda, alpha, gamma, pen, dfmax, eps1, eps2, maxstep1, maxstep2)
        df1 = NULL
        for(k in 1:nlam){
          activeF1 = matrix(fit$betapath[,k],nrow=q)
          df1 = c(df1,median(rowSums(activeF1)))
        }
        df = df1*r1
        loglikelih = (n*q)*log(fit$likhd/(n*q))
        bic <- switch (method,
                       BIC = loglikelih + log(n*q)*df,
                       AIC = loglikelih + 2*df,
                       GCV = fit$likhd*(n*q)/(n*q-df)^2,
                       EBIC = loglikelih + log(n*q)*df + 2*(lgamma(q*p*(p+1)/2+1) 
                                                            - lgamma(df+1) - lgamma(q*p*(p+1)/2-df+1))
        )
      }
      selected = which.min(bic)
      lambda_opt = lambda[selected]
      
      if(isPenColumn){
        fit_opt = EstPenColumn(Y,Z,as.matrix(A),as.matrix(B),as.matrix(C),as.matrix(S),
                               intercept,mu, lambda[1:selected], alpha, gamma, pen, dfmax, eps1, eps2, maxstep1, maxstep2)
        activeF = activeX = fit_opt$betapath[,selected]
      }
      else{
        fit_opt = EstPenSingle(Y,Z,as.matrix(A),as.matrix(B),as.matrix(C),as.matrix(S),
                               intercept,mu,lambda[1:selected], alpha, gamma, pen, dfmax, eps1, eps2, maxstep1, maxstep2)
        activeF = matrix(fit_opt$betapath[,selected],q,p)
        activeX = fit_opt$activeXpath[,selected]
        
      }
      Dnew = fit_opt$Dnew
      mu = fit_opt$mu
    }
    if(method=="CV"&&nlam>1){
      len_cv = ceiling(n/ncv)
      RSS = rep(0,nlam)
      for(jj in 1:ncv){
        cv.id = ((jj-1)*len_cv+1):(jj*len_cv)
        if(jj==ncv) cv.id = ((jj-1)*len_cv+1):n
        Ytrain = Y[-cv.id,]
        Xtrain = X[-cv.id,]
        Ytest = Y[cv.id,]
        Xtest = X[cv.id,]
        Ztrain = bsbasefun(Xtrain,K,degr) 
        Ztest = bsbasefun(Xtest,K,degr)
        if(isPenColumn)
          fit = EstPenColumnCV(Ytrain,Ztrain,Ytest,Ztest,as.matrix(A),as.matrix(B),as.matrix(C),as.matrix(S),
                               intercept,mu,lambda,alpha, gamma, pen, dfmax, eps1,eps2,maxstep1,maxstep2)
        else
          fit = EstPenSingleCV(Ytrain,Ztrain,Ytest,Ztest,as.matrix(A),as.matrix(B),as.matrix(C),as.matrix(S),
                               intercept,mu,lambda,alpha, gamma, pen, dfmax, eps1,eps2,maxstep1,maxstep2)
        RSS = RSS + fit$likhd
      } 
      selected = which.min(RSS)
      lambda_opt = lambda[selected]
      if(isPenColumn){
        fit_opt = EstPenColumn(Y,Z,as.matrix(A),as.matrix(B),as.matrix(C),as.matrix(S),
                               intercept,mu,lambda[1:selected], alpha, gamma, pen, dfmax, eps1, eps2, maxstep1, maxstep2)
        activeF = activeX = fit_opt$betapath[,selected]
      }
      else{
        fit_opt = EstPenSingle(Y,Z,as.matrix(A),as.matrix(B),as.matrix(C),as.matrix(S),
                               intercept,mu,lambda[1:selected], alpha, gamma, pen, dfmax, eps1, eps2, maxstep1, maxstep2)
        activeF = matrix(fit_opt$betapath[,selected],q,p)
        activeX = fit_opt$activeXpath[,selected]
        
      }
      Dnew = fit_opt$Dnew
      mu = fit_opt$mu
     
    }
    if(method=="CV"&&nlam==1){
      if(isPenColumn)
        fit = EstPenColumn(Y,Z,as.matrix(A),as.matrix(B),as.matrix(C),as.matrix(S),
                           intercept,mu,lambda, alpha, gamma, pen, dfmax, eps1, eps2, maxstep1, maxstep2)
      else
        fit = EstPenSingle(Y,Z,as.matrix(A),as.matrix(B),as.matrix(C),as.matrix(S),
                           intercept,mu,lambda, alpha, gamma, pen, dfmax, eps1, eps2, maxstep1, maxstep2)
      selected = 1
      lambda_opt = lambda
      activeX = activeF = fit$betapath
      Dnew = fit$Dnew
      mu = fit$mu
    }

    return(list(Dnew=Dnew,
                betapath=fit$betapath, 
                rss=fit$likhd[selected],
                mu = mu,
                df = fit$df,
                lambda = lambda,
                lambda_opt=lambda_opt,
                selectedID = selected,
                activeF = activeF,
                activeX = activeX,
                Y = Y,
                X = X,
                Z = Z,
                degr = degr,
                K = K
                )
           )
  }