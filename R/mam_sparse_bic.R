
##--------------Estimation with Penalty by BIC----------------------##
mam_sparse_bic <- 
  function(Y,X,method,K_index,r1_index,r2_index,r3_index,pen,isPenColumn,lambda=lambda,A,B,C,S,
           intercept,mu,nlam,degr,lam_min,eps1,maxstep1,eps2,maxstep2,gamma,dfmax,alpha){
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  
  RSS = NULL
  for(K in K_index){
    Z = bsbasefun(X,K,degr)
    for(r3 in r3_index){
      for(r2 in r2_index){
        for(r1 in r1_index){
          if(isPenColumn){
            fit = EstPenColumn(Y,Z,as.matrix(A[,1:r1]),as.matrix(B[1:K,1:r2]),as.matrix(C[,1:r3]),as.matrix(S[1:r3,1:(r1*r2)]),
                               intercept,mu,lambda,alpha,gamma,pen,dfmax,eps1,eps2,maxstep1,maxstep2) 
            df = r1*r2*r3+fit$df*r1+K*r2+q*r3-r1^2-r2^2-r3^2
          }
          else{
            fit = EstPenSingle(Y,Z,as.matrix(A[,1:r1]),as.matrix(B[1:K,1:r2]),as.matrix(C[,1:r3]),as.matrix(S[1:r3,1:(r1*r2)]),
                               intercept,mu,lambda,alpha,gamma,pen,dfmax,eps1,eps2,maxstep1,maxstep2) 
            df1 = NULL
            for(k in 1:nlam){
              activeF1 = matrix(fit$betapath[,k],nrow=q)
              df1 = c(df1,median(rowSums(activeF1)))
            }
            df = r1*r2*r3+df1*r1+K*r2+q*r3-r1^2-r2^2-r3^2
          }
          loglikelih =  -n*q * (log(2*pi) + log(fit$likhd))
          if(method=="BIC"){
            #bic = log(fit$likhd/(n*q)) + log(n*q)*df/(n*q)          
            bic = loglikelih + log(n*q)*df
          }
          if(method=="AIC") bic = loglikelih + 2*df
          if(method=="EBIC"){
            bic = loglikelih + log(n*q)*df
            bic = bic + 2*(lgamma(p+1) - lgamma(df+1) - lgamma(p-df+1))
          }
          if(method=="GCV") bic = loglikelih/(1-df/n)^2
          RSS = cbind(RSS,bic)
        }
      }
    }
  }
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
  if(isPenColumn){
    fit_opt = EstPenColumn(Y,Z,as.matrix(A[,1:r1_opt]),as.matrix(B[1:K_opt,1:r2_opt]),as.matrix(C[,1:r3_opt]),as.matrix(S[1:r3_opt,1:(r1_opt*r2_opt)]),
                           intercept,mu,lambda[1:qj1],alpha, gamma, pen, dfmax, eps1, eps2, maxstep1, maxstep2) 
    activeF = activeX = fit_opt$betapath[,qj1]
  }
  else{
    fit_opt = EstPenSingle(Y,Z,as.matrix(A[,1:r1_opt]),as.matrix(B[1:K_opt,1:r2_opt]),as.matrix(C[,1:r3_opt]),as.matrix(S[1:r3_opt,1:(r1_opt*r2_opt)]),
                           intercept,mu,lambda[1:qj1],alpha, gamma, pen, dfmax, eps1, eps2, maxstep1, maxstep2) 
    activeF = matrix(fit_opt$betapath[,qj1],q,p)
    activeX = fit_opt$activeXpath[,qj1]
  }
  return(list(Dnew=fit_opt$Dnew, 
              rss=fit_opt$likhd[qj1],
              df = fit_opt$df,
              mu = fit_opt$mu,
              activeF = activeF,
              activeX = activeX,
              lambda = lambda,
              selectedID = selected,
              lambda_opt=lambda_opt,
              RSS = RSS,
              rk_opt=c(r1_opt,r2_opt,r3_opt,K_opt),
              Y = Y,
              X = X,
              Z = Z,
              degr = degr,
              K = K
              )
         )
}