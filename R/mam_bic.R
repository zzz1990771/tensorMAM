
##--------------main by BIC without sparsity----------------------##
mam_bic <- function(Y,X,criteria,r1_index,r2_index,r3_index,S,A,B,C,intercept,K,degr,eps,max_step){
  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)
  Ybar = colMeans(Y)
  Y1 = Y - matrix(rep(Ybar,each=n),n)
  Z = bsbasefun(X,K,degr)
  Zbar = colMeans(Z)
  Z = Z - matrix(rep(Zbar,each=n),n)
  RSS = NULL
  for(r3 in r3_index){
    for(r2 in r2_index){
      for(r1 in r1_index){
        fit = Estimation(Y1,Z,as.matrix(A[,1:r1]),as.matrix(B[1:K,1:r2]),as.matrix(C[,1:r3]),
                         as.matrix(S[1:r3,1:(r1*r2)]),eps,max_step)
        df = r1*r2*r3+p*r1+K*r2+q*r3-r1^2-r2^2-r3^2
        loglikelih = (n*q)*log(fit$likhd/(n*q))
        bic <- switch (criteria,
                       BIC = loglikelih + log(n*q)*df,
                       AIC = loglikelih + 2*df,
                       GCV = fit$likhd*(n*q)/(n*q-df)^2,
                       EBIC = loglikelih + log(n*q)*df + 2*(lgamma(q*p*(p+1)/2+1) 
                                                            - lgamma(df+1) - lgamma(q*p*(p+1)/2-df+1))
        )
        RSS = c(RSS,bic)
      }
    }
  }
  selected = which.min(RSS)
  opt = assig(c(length(r1_index),length(r2_index),length(r3_index)))[,selected]
  r1_opt = r1_index[opt[1]]
  r2_opt = r2_index[opt[2]]
  r3_opt = r3_index[opt[3]]
  #---------------- The estimation after selection ---------------------#
  fit = Estimation(Y1,Z,as.matrix(A[,1:r1_opt]),as.matrix(B[,1:r2_opt]),as.matrix(C[,1:r3_opt]),
                   as.matrix(S[1:r3_opt,1:(r1_opt*r2_opt)]),eps,max_step)  
  if(intercept)  mu = Ybar-fit$Dnew%*%Zbar
  else mu = rep(0,q)
  return(list(Dnew=fit$Dnew, 
              rss=fit$likhd,
              mu = mu,
              rk_opt=c(r1_opt,r2_opt,r3_opt),
              selected=selected,
              Y = Y,
              X = X,
              Z = Z,
              degr = degr,
              K = K
              )
         )
}