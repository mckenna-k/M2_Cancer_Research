process <- function(n,l01,l02,l12,censor=NULL){
  if(is.null(censor)){
    t1 <- rexp(n,l01+l02) #temps de passage de 0 à 1 ou 0 à 2
    s1 <- rbinom(n,1,l01/(l01+l02))+1 #détermination de l'état atteint (1 ou 2)
    t2 <- rexp(n,l12) #temps absolu de passage de l'état 1 à 2
    df <- data.frame("DSS"=rep(1,n),"DSS.time"=ifelse(s1==1,t1+t2,t1),
                     "PFI"=rep(1,n),"PFI.time"=t1,
                     "T01"=-(s1-2),"T01.time"=t1,
                     "T02"=s1-1,"T02.time"=ifelse(s1==2,t1,t1+t2),
                     "T12"=-(s1-2),"T12.time"=ifelse(s1==1,t2,0))
  }
  else{
    t1 <- rexp(n,l01+l02+censor) #temps de passage de 0 à 1 ou 0 à 2
    s1 <- apply(rmultinom(n, size = 1, prob = c(censor/(l01+l02+censor),l01/(l01+l02+censor),l02/(l01+l02+censor)))==1,2,which)-1 #détermination de l'état atteint (1 ou 2)
    t2 <- rexp(n,l12+censor)
    s2 <- rbinom(n,1,censor/(l01+censor))
    df <- data.frame("DSS"=ifelse(s1 > 0,1,0),"DSS.time"=ifelse(s1==1,t1+t2,t1),
                     "PFI"=ifelse(s1 > 0,1,0),"PFI.time"=t1,
                     "T01"=ifelse(s1==1,1,0),"T01.time"=t1,
                     "T02"=ifelse(s1-1 > 0,1,0),"T02.time"=ifelse(s1==2 | s1==0,t1,t1+t2),
                     "T12"=ifelse(s1==1 & s2==1,1,0),"T12.time"=ifelse(s1==1,t2,0))
  }
  return(df)
}

df_new <- process(1000,0.3,0.3,0.3,0.9) #fonction générant la data frame avec une intensité de censure grande pour en avoir beaucoup
covar <- sample(c(0,1),11600,replace=TRUE)
covar_sim <- sample(c(0,1),1000,replace=TRUE)