
##--------------Estimation without Penalty----------------------##
mam_cv <- function(Y,X,ncv,r1_index,r2_index,r3_index,S,A,B,C,intercept,K,degr,eps,max_step){
  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)
  Ybar = colMeans(Y)
  Y1 = Y - matrix(rep(Ybar,each=n),n)
  Z = bsbasefun(X,K,degr) 
  Zbar = colMeans(Z)
  Z = Z - matrix(rep(Zbar,each=n),n)
  
  len_cv = floor(n/ncv)  
  RSS = rep(0,length(r1_index)*length(r2_index)*length(r3_index))
  for(jj in 1:ncv){ # start CV
    cv.id = ((jj-1)*len_cv+1):(jj*len_cv)
    if(jj==ncv) cv.id = ((jj-1)*len_cv+1):n
    Ytrain = Y1[-cv.id,]
    Ztrain = Z[-cv.id,]
    Ytest = Y1[cv.id,]
    Ztest = Z[cv.id,]
    
    RSS0 = NULL
    for(r3 in r3_index){
      for(r2 in r2_index){
        for(r1 in r1_index){
          fit = Estimation(Ytrain,Ztrain,as.matrix(A[,1:r1]),as.matrix(B[1:K,1:r2]),as.matrix(C[,1:r3]),
                           as.matrix(S[1:r3,1:(r1*r2)]),eps,max_step)
          Dnew  = fit$Dnew
          RSS0 = c(RSS0,sum((Ytest - Ztest%*%t(Dnew))^2))
        }
      }
    }
    RSS = RSS + RSS0
  } # end of CV
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
