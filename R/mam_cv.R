
##--------------Estimation without Penalty----------------------##
mam_cv <- function(Y,X,ncv,K_index,r1_index,r2_index,r3_index,S,A,B,C,degr,eps,max_step){
  n <- dim(Y)[1]
  q <- dim(Y)[2]
  p <- dim(X)[2]
  len_cv = ceiling(n/ncv)  
  RSS = rep(0,length(r1_index)*length(r2_index)*length(r3_index)*length(K_index))
  for(jj in 1:ncv){ # start CV
    cv.id = ((jj-1)*len_cv+1):(jj*len_cv)
    if(jj==ncv) cv.id = ((jj-1)*len_cv+1):n
    Ytrain = Y[-cv.id,]
    Xtrain = X[-cv.id,]
    Ytest = Y[cv.id,]
    Xtest = X[cv.id,]
    
    RSS0 = NULL
    for(K in K_index){
      Ztrain = bsbasefun(Xtrain,K,degr)
      Ztest = bsbasefun(Xtest,K,degr)
      for(r3 in r3_index){
        for(r2 in r2_index){
          for(r1 in r1_index){
            fit = Estimation(Ytrain,Ztrain,as.matrix(A[,1:r1]),as.matrix(B[1:K,1:r2]),as.matrix(C[,1:r3]),as.matrix(S[1:r3,1:(r1*r2)]),eps,max_step)
            Dnew  = fit$Dnew
            RSS0 = c(RSS0,sum((Ytest - Ztest%*%t(Dnew))^2))
          }
        }
      }
    }
    RSS = RSS + RSS0
  } # end of CV
  selected = which.min(RSS)
  opt = assig(c(length(r1_index),length(r2_index),length(r3_index),length(K_index)))[,selected]
  r1_opt = r1_index[opt[1]]
  r2_opt = r2_index[opt[2]]
  r3_opt = r3_index[opt[3]]
  K_opt = K_index[opt[4]]
  #---------------- The estimation after selection ---------------------#
  Z = bsbasefun(X,K_opt,degr)
  fit = Estimation(Y,Z,as.matrix(A[,1:r1_opt]),as.matrix(B[1:K_opt,1:r2_opt]),as.matrix(C[,1:r3_opt]),as.matrix(S[1:r3_opt,1:(r1_opt*r2_opt)]),eps,max_step)
  return(list(Dnew=fit$Dnew, 
              rss=fit$likhd,
              rk_opt=c(r1_opt,r2_opt,r3_opt,K_opt),
              selected=selected,
              Y = Y,
              X = X,
              Z = Z,
              degr = degr,
              K = K
              )
         )
}
