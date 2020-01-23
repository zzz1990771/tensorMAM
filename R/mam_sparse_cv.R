
##--------------Estimation with Penalty by CV----------------------##
mam_sparse_cv <- 
  function(Y,X,ncv,K_index,r1_index,r2_index,r3_index,pen,isPenColumn,lambda=lambda,A,B,C,S,
           intercept,nlam,degr,lam_min,eps1,maxstep1,eps2,maxstep2,gamma,dfmax,alpha){
    n <- dim(Y)[1]
    q <- dim(Y)[2]
    p <- dim(X)[2]
    Ybar = colMeans(Y)
    Y1 = Y - matrix(rep(Ybar,each=n),n)
    
    len_cv = ceiling(n/ncv)
    RSS = matrix(0,nlam,length(r1_index)*length(r2_index)*length(r3_index)*length(K_index))
    for(jj in 1:ncv){ # start CV
      cv.id = ((jj-1)*len_cv+1):(jj*len_cv)
      if(jj==ncv) cv.id = ((jj-1)*len_cv+1):n
      Ytrain = Y1[-cv.id,]
      Xtrain = X[-cv.id,]
      Ytest = Y1[cv.id,]
      Xtest = X[cv.id,]
      RSS0 = NULL
      for(K in K_index){
        Ztrain = bsbasefun(Xtrain,K,degr) 
        Ztest = bsbasefun(Xtest,K,degr)
        for(r3 in r3_index){
          for(r2 in r2_index){
            for(r1 in r1_index){
              if(isPenColumn)
                fit = EstPenColumnCV(Ytrain,Ztrain,Ytest,Ztest,as.matrix(A[,1:r1]),as.matrix(B[1:K,1:r2]),as.matrix(C[,1:r3]),as.matrix(S[1:r3,1:(r1*r2)]),
                                     lambda,alpha, gamma, pen, dfmax, eps1,eps2,maxstep1,maxstep2) 
              else
                fit = EstPenSingleCV(Ytrain,Ztrain,Ytest,Ztest,as.matrix(A[,1:r1]),as.matrix(B[1:K,1:r2]),as.matrix(C[,1:r3]),as.matrix(S[1:r3,1:(r1*r2)]),
                                     lambda,alpha, gamma, pen, dfmax, eps1,eps2,maxstep1,maxstep2) 
              RSS0 = cbind(RSS0,fit$likhd)
            }
          }
        }
      }
      RSS = RSS + RSS0
    } # end of CV
    selected = which.min(RSS)
    qj = ceiling(selected/nlam)
    qj1 = selected%%nlam
    if(qj1==0) qj1=nlam
    
    lambda_opt = lambda[qj1]
    opt = assig(c(length(r1_index),length(r2_index),length(r3_index),length(K_index)))[,qj]
    r1_opt = r1_index[opt[1]]
    r2_opt = r2_index[opt[2]]
    r3_opt = r3_index[opt[3]]
    K_opt = K_index[opt[4]]
    
    #---------------- The estimation after selection ---------------------#
    Z = bsbasefun(X,K_opt,degr)
    Zbar = colMeans(Z)
    Z = Z - matrix(rep(Zbar,each=n),n)
    if(isPenColumn){
      fit_opt = EstPenColumn(Y1,Z,as.matrix(A[,1:r1_opt]),as.matrix(B[1:K_opt,1:r2_opt]),as.matrix(C[,1:r3_opt]),as.matrix(S[1:r3_opt,1:(r1_opt*r2_opt)]),
                             lambda[1:qj1],alpha, gamma, pen, dfmax,eps1, eps2, maxstep1, maxstep2) 
      activeF = activeX = fit_opt$betapath[,qj1]
    }
    else{ 
      fit_opt = EstPenSingle(Y1,Z,as.matrix(A[,1:r1_opt]),as.matrix(B[1:K_opt,1:r2_opt]),as.matrix(C[,1:r3_opt]),as.matrix(S[1:r3_opt,1:(r1_opt*r2_opt)]),
                             lambda[1:qj1],alpha, gamma, pen, dfmax,eps1, eps2, maxstep1, maxstep2) 
      activeF = matrix(fit_opt$betapath[,qj1],q,p)
      activeX = fit_opt$activeXpath[,qj1]
    }
    if(intercept)  mu = Ybar-Dnew%*%Zbar
    else mu = rep(0,q)
    return(list(Dnew=fit_opt$Dnew, 
                rss=fit_opt$likhd,
                df = fit_opt$df,
                mu = mu,
                activeF = activeF,
                activeX = activeX,
                lambda = lambda,
                selectedID = selected,
                RSS = RSS,
                lambda_opt=lambda_opt,
                rk_opt=c(r1_opt,r2_opt,r3_opt,K_opt),
                Y = Y,
                X = X,
                Z = Z,
                degr = degr,
                K = K
    )
    )
  }