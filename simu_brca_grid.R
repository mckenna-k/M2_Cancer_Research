source("common.R")

library(msm)

# 2: rapport d'individus en stage III-IV parmi ceux non censurés dans la transi 0 -> 1
# 3: rapport d'individus en stage III-IV parmi ceux non censurés dans la transi 0 -> 2
# 4: nombre d'individus en stage III-IV parmi ceux non censurés dans la transi 1 -> 2
# 5: nombre d'individus non censurés dans la transi 1 -> 2
# 6: nombre d'individus en stage III-IV dans le cancer
covar_creation <- function(df_new,p1,p2,p3num,p3denom,total,seed=1){
  set.seed(seed)
  covar <- rep(0,nrow(df_new))
  covar[df_new$T01==1] <- sample(0:1,sum(df_new$T01==1),replace = TRUE,prob = c(1-p1,p1))
  covar[df_new$T02==1] <- sample(0:1,sum(df_new$T02==1),replace = TRUE,prob = c(1-p2,p2))
  deja_pris <- sum(df_new$T01==1 & covar==1)
  if(deja_pris < 25) covar[df_new$T12==1 & df_new$T01==0] <- sample(0:1,size=sum(df_new$T12==1)-deja_pris,replace = TRUE,prob = c(1-(p3num-deja_pris)/(p3denom-deja_pris),(p3num-deja_pris)/(p3denom-deja_pris)))
  non_pris <- sum(df_new$T01==0 & df_new$T02==0 & df_new$T12==0)
  if(total-sum(covar) > 0) covar[df_new$T01==0 & df_new$T02==0 & df_new$T12==0] <- sample(0:1,non_pris,replace = TRUE,prob=c(1-(total-sum(covar))/non_pris,(total-sum(covar))/non_pris))
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

params = expand.grid(
  i01 = seq(0.1,0.9,0.4),
  i02 = seq(0.1,0.9,0.4),
  i12 = seq(0.1,0.9,0.4),
  seeds = 1:3
)

# val01 = .6; val02 = .2; val12 = .3; # BRCA
# val01 = .03; val02 = .01; val12 = .2; # BRCA crud init
# val01 = .7; val02 = .1; val12 = .5; # OV

head(params)

res = apply(params, 1, function(p) {
  i01 = p[[1]]
  i02 = p[[2]]
  i12 = p[[3]]
  seed =p[[4]]
  print(paste("******", i01,i02,i12, seed))
  df_new <- 
  df_msm <- df_msm_creation(df_new)
  covar <- covar_creation(df_new,42/98,16/33,25,46,257,seed=seed)
  df_msm[["covar"]] <- covar[df_msm$indiv]
  msm_simu <- msm(state~time, subject=indiv,data=df_msm, gen.inits=TRUE, qmatrix= Qmat, deathexact = 3, covariates = ~covar)
  hazards <- msm::hazard.msm(msm_simu)[[1]]
  if(dim(hazards)[2] == 1){
    deltahr01 = 10000
    deltahr02 = 10000
    deltahr03 = 10000
  } else{
    deltahr01 =  hazards[1,3]-hazards[1,2]
    deltahr02 =  hazards[2,3]-hazards[2,2]
    deltahr03 =  hazards[3,3]-hazards[3,2]
  }
  ret = c(i01=i01, i02=i02, i12=i12, seed=seed, deltahr01=deltahr01, deltahr02=deltahr02, deltahr03=deltahr03)
  return(ret)
})




d = data.frame(t(res))
head(d)

max_deltahr = 100
d[d$deltahr01 > max_deltahr,]$deltahr01 = max_deltahr
d[d$deltahr02 > max_deltahr,]$deltahr02 = max_deltahr
d[d$deltahr03 > max_deltahr,]$deltahr03 = max_deltahr

d$i01 = as.factor(d$i01)
d$i02 = as.factor(d$i02)
d$i12 = as.factor(d$i12)


layout(matrix(1:3,1), respect=TRUE)
boxplot(deltahr01~i01*i02*i12, d, las=2, xlab="")
boxplot(deltahr02~i01*i02*i12, d, las=2, xlab="")
boxplot(deltahr03~i01*i02*i12, d, las=2, xlab="")
par(mar=c(5.1, 4.1, 4.1, 2.1))





stop()



hist(d$deltahr01)
plot(density(unlist(d$deltahr01)))
 
m = lm(deltahr01~i01*i02*i12 , d)
anova(m)

head(d)

d$i01 = as.factor(d$i01)
d$i02 = as.factor(d$i02)
d$i12 = as.factor(d$i12)

m = lm(deltahr01~i01*i02*i12 , d)
anova(m)

par(mar=c(10, 4.1, 4.1, 2.1))
boxplot(deltahr01~i01*i02*i12, d, las=2, xlab="")
par(mar=c(5.1, 4.1, 4.1, 2.1))


d[d$deltahr01 > 100,]$deltahr01 = 100
d[d$deltahr02 > 100,]$deltahr02 = 100
d[d$deltahr03 > 100,]$deltahr03 = 100

m = lm(deltahr01~i01*i02*i12 , d)
anova(m)

# BRCA 0.6, 0.2, 0.3
# OV 0.7, 0.1, 0.5
par(mar=c(10, 4.1, 4.1, 2.1))












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

# obs: diminuer le rapport i01/i02 augmente instabilité ?