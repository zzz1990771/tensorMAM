
##--------------Estimation with Penalty by BIC----------------------##
mam_sparse_bic <- 
  function(Y,X,K_index,r1_index,r2_index,r3_index,pen,lambda=lambda,A,B,C,S,
           nlam,degr,lam_min,eps1,maxstep1,eps2,maxstep2,gamma,dfmax,alpha,setlam){
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  
  RSS = NULL
  activeA = NULL
  for(K in K_index){
    Z = bsbasefun(X,K,degr)
    for(r3 in r3_index){
      for(r2 in r2_index){
        for(r1 in r1_index){
          fit = Estimation_penalty(Y,Z,as.matrix(A[,1:r1]),as.matrix(B[1:K,1:r2]),as.matrix(C[,1:r3]),as.matrix(S[1:r3,1:(r1*r2)]),
                                   lambda,alpha,gamma,pen,dfmax,eps1,eps2,maxstep1,maxstep2,0,setlam) 
          df = r1*r2*r3+fit$df*r1+K*r2+q*r3-r1^2-r2^2-r3^2
          bic = 2*log(fit$likhd) + log(n)*df/n
          RSS = cbind(RSS,bic)
          activeA = cbind(activeA,fit$betapath)
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
  fit_opt = Estimation_penalty(Y,Z,as.matrix(A[,1:r1_opt]),as.matrix(B[1:K_opt,1:r2_opt]),as.matrix(C[,1:r3_opt]),as.matrix(S[1:r3_opt,1:(r1_opt*r2_opt)]),
                           lambda[1:qj1],alpha, gamma, pen, dfmax, eps1, eps2, maxstep1, maxstep2,0,setlam) 
  return(list(Dnew=fit_opt$Dnew, 
              rss=fit_opt$likhd[qj1],
              df = fit_opt$df,
              activeA = fit_opt$betapath[,qj1],
              lambda = lambda,
              selectedID = selected,
              lambda_opt=lambda_opt,
              RSS = RSS,
              rk_opt=c(r1_opt,r2_opt,r3_opt,K_opt),
              Y = Y,
              X = X,
              Z = Z
              )
         )
}