library(msm)

process <- function(n,l01,l02,l12,censor1=NULL,censor2=NULL,seed=1){
  set.seed(seed)
  if(is.null(censor1) | is.null(censor2)){
    t1 <- rexp(n,l01+l02) #temps de passage de 0 à 1 ou 0 à 2
    s1 <- rbinom(n,1,l02/(l01+l02))+1 #détermination de l'état atteint (1 ou 2)
    t2 <- rexp(n,l12) #temps de passage de l'état 1 à 2 à partir du temps d'arrivée à 1
    df <- data.frame("DSS"=rep(1,n),"DSS.time"=ifelse(s1==1,t1+t2,t1),
                     "PFI"=rep(1,n),"PFI.time"=t1,
                     "T01"=-(s1-2),"T01.time"=t1,
                     "T02"=s1-1,"T02.time"=ifelse(s1==2,t1,t1+t2),
                     "T12"=-(s1-2),"T12.time"=ifelse(s1==1,t2,0))
  }
  else{
    t1 <- rexp(n,l01+l02+censor1) #temps de passage de 0 à 1 ou 0 à 2
    s1 <- apply(rmultinom(n, size = 1, prob = c(censor1/(l01+l02+censor1),l01/(l01+l02+censor1),l02/(l01+l02+censor1)))==1,2,which)-1 #détermination de l'état atteint (1 ou 2)
    t2 <- rexp(n,l12+censor2)
    s2 <- rbinom(n,1,l12/(l12+censor2))
    df <- data.frame("DSS"=ifelse(s1 > 0,1,0),"DSS.time"=ifelse(s1==1,t1+t2,t1),
                     "PFI"=ifelse(s1 > 0,1,0),"PFI.time"=t1,
                     "T01"=ifelse(s1==1,1,0),"T01.time"=t1,
                     "T02"=ifelse(s1-1 > 0,1,0),"T02.time"=ifelse(s1==2 | s1==0,t1,t1+t2),
                     "T12"=ifelse(s1==1 & s2==1,1,0),"T12.time"=ifelse(s1==1,t2,0))
  }
  return(df)
}

covar_creation <- function(df_new,p1,p2,p3num,p3denom,seed=1){
  set.seed(seed)
  covar <- rep(0,nrow(df_new))
  covar[df_new$T01==1] <- sample(0:1,sum(df_new$T01==1),replace = TRUE,prob = c(1-p1,p1))
  covar[df_new$T02==1] <- sample(0:1,sum(df_new$T02==1),replace = TRUE,prob = c(1-p2,p2))
  deja_pris <- sum(df_new$T01==1 & covar==1)
  if(deja_pris < sum(df_new$T12==1) & deja_pris < 25) covar[df_new$T12==1 & df_new$T01==0] <- sample(0:1,size=sum(df_new$T12==1)-deja_pris,replace = TRUE,prob = c(1-(p3num-deja_pris)/(p3denom-deja_pris),(p3num-deja_pris)/(p3denom-deja_pris)))
  return(covar)
}

df_msm_creation <- function(df_new){
  indiv_msm <- c((1:nrow(df_new)),
                 (1:nrow(df_new))[(df_new$T01.time < df_new$T02.time) & df_new$T01==1],
                 (1:nrow(df_new))[df_new$T12.time > 0 & df_new$T12==1],
                 (1:nrow(df_new))[df_new$T12.time > 0 & df_new$T12==0],
                 (1:nrow(df_new))[df_new$T02 == 1],
                 (1:nrow(df_new))[df_new$T02 == 0 & df_new$T01==0])
  state_msm <- c(rep(1,nrow(df_new)),
                 rep(2,sum((df_new$T01.time < df_new$T02.time) & df_new$T01==1)),
                 rep(3,sum(df_new$T12.time > 0 & df_new$T12==1)),
                 rep(2,sum(df_new$T12.time > 0 & df_new$T12==0)),
                 rep(3,sum(df_new$T02 == 1)),
                 rep(1,sum(df_new$T02 == 0 & df_new$T01==0)))
  time_msm <- c(rep(0,nrow(df_new)),
                df_new$T01.time[(df_new$T01.time < df_new$T02.time) & df_new$T01==1],
                (df_new$T12.time + df_new$T01.time)[df_new$T12.time > 0 & df_new$T12==1],
                (df_new$T12.time + df_new$T01.time)[df_new$T12.time > 0 & df_new$T12==0],
                df_new$T02.time[df_new$T02 == 1],
                df_new$T02.time[df_new$T02 == 0 & df_new$T01==0])
  df_msm <- data.frame("indiv"=indiv_msm,"state"=state_msm,"time"=time_msm)
  df_msm <- df_msm[order(df_msm$indiv,df_msm$time),]
  rownames(df_msm) <- 1:nrow(df_msm)
  return(df_msm)
}

n <- 1000
Qmat <- matrix(c(1,1,1,0,1,1,0,0,0),ncol=3,byrow = TRUE)
val01 <- seq(0.1,0.9,0.2)
val02 <- seq(0.1,0.9,0.2)
val12 <- seq(0.1,0.9,0.2)

res01 <-rep(NA,length(val01))
res02 <-rep(NA,length(val02))
res12 <-rep(NA,length(val12))

#ajouter paramètre set.seed()
#grid.search()
#10% de censure
#outcome : taille de l'intervalle de confiance
#grille tous les seq(0.1,0.9,0.2)  125 points 
#réplicats 3 fois
#fonction covar seed
for(seed in 1:1){
  for(i in 1:length(val01)){
    for(j in 1:length(val02)){
      for(k in 1:length(val12)){
        df_new <- process(n=n,val01[[i]],val02[[j]],val12[[k]],0.001,0.35,seed)
        df_msm <- df_msm_creation(df_new)
        covar <- covar_creation(df_new,42/98,16/33,25,46,seed)
        df_msm[["covar"]] <- covar[df_msm$indiv]
        msm_simu <- msm(state~time, subject=indiv,data = df_msm, gen.inits = TRUE,qmatrix= Qmat, deathexact = 3,covariates = ~covar)
        hazards <- hazard.msm(msm_simu)[[1]]
        if(dim(hazards)[2] == 1){
          res01[k+j*length(val02)+i*length(val02)*length(val01)] <- 10000
          res02[k+j*length(val02)+i*length(val02)*length(val01)] <- 10000
          res12[k+j*length(val02)+i*length(val02)*length(val01)] <- 10000
        }
        else{
          res01[k+j*length(val02)+i*length(val02)*length(val01)] <- hazards[1,3]-hazards[1,2]
          res02[k+j*length(val02)+i*length(val02)*length(val01)] <- hazards[2,3]-hazards[2,2]
          res12[k+j*length(val02)+i*length(val02)*length(val01)] <- hazards[3,3]-hazards[3,2]
        }
      }  
    }
  }
}

mat_res01 <- matrix(res01,byrow = TRUE,ncol=3)[order(x),]
mat_res02 <- matrix(res02,byrow = TRUE,ncol=3)[order(x),]
mat_res12 <- matrix(res12,byrow = TRUE,ncol=3)[order(x),]
x <- x[order(x)]

par(mfrow=c(1,3),mar=rep(2,4))
plot(x,mat_res01[,1],type="l",ylim=c(0,10))
points(x,mat_res01[,2],type="l",col=2,lty=2)
points(x,mat_res01[,3],type="l",col=2,lty=2)
abline(h=1,lty=2,col=3)

plot(x,mat_res02[,1],type="l",ylim=c(0,10))
points(x,mat_res02[,2],type="l",col=2,lty=2)
points(x,mat_res02[,3],type="l",col=2,lty=2)
abline(h=1,lty=2,col=3)

plot(x,mat_res12[,1],type="l",ylim=c(0,10))
points(x,mat_res12[,2],type="l",col=2,lty=2)
points(x,mat_res12[,3],type="l",col=2,lty=2)
abline(h=1,lty=2,col=3)

#pour BRCA T01 <- 98
# T02 <- 33
# T12 <- 46
# T01 III- IV <- 42
# T02 III- IV <- 16
# T12 III- IV <- 25

# df_new <- process(n=n,0.6,0.2,0.3,0.01,0.35,1)
# sum(df_new$T01==1)
# sum(df_new$T02==1)
# sum(df_new$T12==1)

# obs: diminuer le rapport l01/l02 augmente instabilité ?