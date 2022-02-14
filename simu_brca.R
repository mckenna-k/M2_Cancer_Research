library(msm)
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

process <- function(n, i01, i02, i12, i00=NULL, seed=1){
  set.seed(seed)
  if(is.null(i00)){
    t1 <- rexp(n,i01+i02) #temps de passage de 0 à 1 ou 0 à 2
    s1 <- rbinom(n,1,i02/(i01+i02))+1 #détermination de l'état atteint (1 ou 2)
    t2 <- rexp(n,i12) #temps de passage de l'état 1 à 2 à partir du temps d'arrivée à 1
    df <- data.frame("DSS"=rep(1,n),"DSS.time"=ifelse(s1==1,t1+t2,t1),
                     "PFI"=rep(1,n),"PFI.time"=t1,
                     "T01"=-(s1-2),"T01.time"=t1,
                     "T02"=s1-1,"T02.time"=ifelse(s1==2,t1,t1+t2),
                     "T12"=-(s1-2),"T12.time"=ifelse(s1==1,t2,0))
  }  else {
    tcensure <- rexp(n,i00) #temps de censure pour chaque individu
    tmin0102 <- rexp(n,i01+i02) #min temps de passage de 0 à 1 et de 0 à 2
    t1 <- apply(matrix(c(tcensure,tmin0102),byrow = FALSE,ncol=2),1,min) #min temps de passage de 0 à 1 et de 0 à 2 et du temps de censure
    s1 <- ifelse(tmin0102 < tcensure,rbinom(n, 1, i02/(i01+i02))+1,0) #détermination de l'état atteint (1 ou 2 ou censure)
    tmin12 <- rexp(n,i12)
    t2 <- apply(matrix(c(tmin12,tcensure-t1),byrow = FALSE,ncol=2),1,min) #min temps de passage de 1 à 2 et du temps de censure
    s2 <- ifelse(tmin12 < tcensure,1,0) #détermination de la censure de la transition 1 -> 2
    df <- data.frame("DSS"=ifelse(s1 == 2 | (s1 == 1 & s2 == 1),1,0),"DSS.time"=ifelse(s1==1,t1+t2,t1),
                     "PFI"=ifelse(s1 > 0,1,0),"PFI.time"=t1,
                     "T01"=ifelse(s1==1,1,0),"T01.time"=t1,
                     "T02"=ifelse(s1-1 > 0,1,0),"T02.time"=ifelse(s1==2 | s1==0,t1,t1+t2),
                     "T12"=ifelse(s1==1 & s2==1,1,0),"T12.time"=ifelse(s1==1,t2,0))
  }
  return(df)
}




n <- 1000
Qmat <- matrix(c(1,1,1,0,1,1,0,0,0),ncol=3,byrow = TRUE)
val <- seq(10,0.01,length.out=100)
# val <- seq(1,0.01,length.out=100)
contingent01 <- rep(NA,length(val))
contingent02 <- rep(NA,length(val))
contingent12 <- rep(NA,length(val))
res01 <-rep(NA,3*length(val))
res02 <-rep(NA,3*length(val))
res12 <-rep(NA,3*length(val))

covar <- sample(0:1,1000,replace=TRUE,prob=c(0.8,0.2))

for(i in 1:length(val)){
  print(i)
  df_new <- process(n=n,0.6,0.2,0.3,val[[i]],i)
  df_msm <- df_msm_creation(df_new)
  df_msm[["covar"]] <- covar[df_msm$indiv]
  msm_simu <- msm(state~time, subject=indiv,data = df_msm, gen.inits = TRUE,qmatrix= Qmat, deathexact = 3,covariates = ~covar)
  hazards <- hazard.msm(msm_simu)[[1]]
  if(dim(hazards)[2] == 1){
    res01[((i-1)*3+1):(i*3)] <- c(hazards[1,],NA,NA)
    res02[((i-1)*3+1):(i*3)] <- c(hazards[2,],NA,NA)
    res12[((i-1)*3+1):(i*3)] <- c(hazards[3,],NA,NA)
  }
  else{
    res01[((i-1)*3+1):(i*3)] <- hazards[1,]
    res02[((i-1)*3+1):(i*3)] <- hazards[2,]
    res12[((i-1)*3+1):(i*3)] <- hazards[3,]
  }
  contingent01[[i]] <- sum(df_new$T01==1)
  contingent02[[i]] <- sum(df_new$T02==1)
  contingent12[[i]] <- sum(df_new$T12==1)
}

x <- contingent01+contingent02
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

plot(val,contingent01/contingent02,type="l",ylim=c(0,5))
points(val,contingent01/contingent12,type="l",col=2)
points(val,contingent12/contingent02,type="l",col=3)
abline(h=98/33,lty=2,col=1)
abline(h=98/46,lty=2,col=2)
abline(h=46/33,lty=2,col=3)
legend("bottomleft",legend=c("3","2.1","1.4"),lty=2,col=1:3)

#pour BRCA T01 <- 98
# T02 <- 33
# T12 <- 46
# T01 III- IV <- 42
# T02 III- IV <- 16
# T12 III- IV <- 25

df_new <- process(n=n,0.6,0.2,0.3,6)
sum(df_new$T01==1)
sum(df_new$T02==1)
sum(df_new$T12==1)

# obs: diminuer le rapport l01/l02 augmente instabilité ?