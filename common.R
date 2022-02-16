if (!exists("mmsm")) {mmsm = memoise::memoise(function(formula, subject, ...) msm::msm(formula, subject, ...))}


# process.orig <- function(n, i01, i02, i12, censor_param=NULL, seed=1){
#   set.seed(seed)
#   if(is.null(censor_param)){
#     t1 <- rexp(n,i01+i02) #temps de passage de 0 à 1 ou 0 à 2
#     s1 <- rbinom(n,1,i02/(i01+i02))+1 #détermination de l'état atteint (1 ou 2)
#     t2 <- rexp(n,i12) #temps de passage de l'état 1 à 2 à partir du temps d'arrivée à 1
#     dfsimu <- data.frame("DSS"=rep(1,n),"DSS.time"=ifelse(s1==1,t1+t2,t1),
#                      "PFI"=rep(1,n),"PFI.time"=t1,
#                      "T01"=-(s1-2),"T01.time"=t1,
#                      "T02"=s1-1,"T02.time"=ifelse(s1==2,t1,t1+t2),
#                      "T12"=-(s1-2),"T12.time"=ifelse(s1==1,t2,0))
#   }  else {
#     tcensure <- rexp(n,censor_param) #temps de censure pour chaque individu
#     tmin0102 <- rexp(n,i01+i02) #min temps de passage de 0 à 1 et de 0 à 2
#     t1 <- apply(matrix(c(tcensure,tmin0102),byrow = FALSE,ncol=2),1,min) #min temps de passage de 0 à 1 et de 0 à 2 et du temps de censure
#     s1 <- ifelse(tmin0102 < tcensure,rbinom(n, 1, i02/(i01+i02))+1,0) #détermination de l'état atteint (1 ou 2 ou censure)
#     tmin12 <- rexp(n,i12)
#     t2 <- apply(matrix(c(tmin12,tcensure-t1),byrow = FALSE,ncol=2),1,min) #min temps de passage de 1 à 2 et du temps de censure
#     s2 <- ifelse(tmin12 < tcensure,1,0) #détermination de la censure de la transition 1 -> 2
#     dfsimu <- data.frame("DSS"=ifelse(s1 == 2 | (s1 == 1 & s2 == 1),1,0),"DSS.time"=ifelse(s1==1,t1+t2,t1),
#                      "PFI"=ifelse(s1 > 0,1,0),"PFI.time"=t1,
#                      "T01"=ifelse(s1==1,1,0),"T01.time"=t1,
#                      "T02"=ifelse(s1-1 > 0,1,0),"T02.time"=ifelse(s1==2 | s1==0,t1,t1+t2),
#                      "T12"=ifelse(s1==1 & s2==1,1,0),"T12.time"=ifelse(s1==1,t2,0))
#   }
#   return(dfsimu)
# }

# process(p[1], p[2], p[3], p[4])
#
# n  =p[1]
# i01=p[2]
# i02=p[3]
# i12=p[4]
process <- function(n, i01, i02, i12, censor_param=NULL, seed=1, version="T02timePFI"){
  set.seed(seed)
  t1 <- rexp(n,i01+i02) #temps de passage de 0 à 1 ou 0 à 2
  s1 <- rbinom(n,1,i02/(i01+i02))+1 #détermination de l'état atteint (1 ou 2)
  t2 <- rexp(n,i12) #temps de passage de l'état 1 à 2 à partir du temps d'arrivée à 1
  if(is.null(censor_param)){
    # which tiome for T02.time
    if (version=="T02timeDSS") {
      T02.time = ifelse(s1==2,t1,t1+t2)
    } else {
      T02.time = t1
    }    
    dfsimu <- data.frame("DSS"=rep(1,n),"DSS.time"=ifelse(s1==1,t1+t2,t1),
                     "PFI"=rep(1,n),"PFI.time"=t1,
                     "T01"=-(s1-2),"T01.time"=t1,
                     "T02"=s1-1,"T02.time"=T02.time,
                     "T12"=-(s1-2),"T12.time"=ifelse(s1==1,t2,0))
  }  else {  
    set.seed(seed*2+1)
    tcensure <- rexp(n,censor_param) #temps de censure pour chaque individu
    print(paste0("nb_censored_DSS = ", sum(((t1+t2 > tcensure) & (s1==1)) | ((t1 > tcensure) & (s1==2)) )))
    print(paste0("nb_censored_PFI = ", sum(t1 > tcensure)))
    tmin0102 = t1
    t1 <- apply(matrix(c(tcensure,tmin0102),byrow = FALSE,ncol=2),1,min) #min temps de passage de 0 à 1 et de 0 à 2 et du temps de censure
    s1 <- ifelse(tmin0102 < tcensure,s1,0) #détermination de l'état atteint (1 ou 2 ou censure)
    tmin12 <- t2
    t2 <- apply(matrix(c(tcensure-t1,tmin12),byrow = FALSE,ncol=2),1,min) #min temps de passage de 1 à 2 et du temps de censure
    s2 <- ifelse(tmin12 < tcensure-t1,1,0) #détermination de la censure de la transition 1 -> 2
    if (version=="T02timeDSS") {
      T02.time = ifelse(s1==2 | s1==0,t1,t1+t2)
    } else {
      T02.time = t1
    }    
    dfsimu <- data.frame("DSS"=ifelse(s1 == 2 | (s1 == 1 & s2 == 1),1,0),"DSS.time"=ifelse(s1==1,t1+t2,t1),
                     "PFI"=ifelse(s1 > 0,1,0),"PFI.time"=t1,
                     "T01"=ifelse(s1==1,1,0),"T01.time"=t1,
                     "T02"=ifelse(s1-1 > 0,1,0),"T02.time"=T02.time,
                     "T12"=ifelse(s1==1 & s2==1,1,0),"T12.time"=ifelse(s1==1,t2,0))
  }
  return(dfsimu)
}



add_censor <- function(dfsimu, censor_param, seed=1){
  set.seed(seed*2+1)
  n = nrow(dfsimu)
  tcensure <- rexp(n,censor_param) #temps de censure pour chaque individu
  head(dfsimu)
  # plot(density(df$PFI.time))
  # lines(density(df$DSS.time), col=2)
  # lines(density(tcensure), col=4)
  
  print(paste0("nb_censored_DSS = ", sum(dfsimu$DSS.time > tcensure)))
  print(paste0("nb_censored_PFI = ", sum(dfsimu$PFI.time > tcensure)))

  dfsimu$DSS      = ifelse(dfsimu$DSS.time >= tcensure, 0       , 1)
  dfsimu$DSS.time = ifelse(dfsimu$DSS.time >= tcensure, tcensure, dfsimu$DSS.time)
  dfsimu$PFI      = ifelse(dfsimu$PFI.time >= tcensure, 0       , 1)
  dfsimu$PFI.time = ifelse(dfsimu$PFI.time >= tcensure, tcensure, dfsimu$PFI.time)
  # print(dfsimu$DSS)
  dfsimu$T01 = dfsimu$T01.time = dfsimu$T02 = dfsimu$T02.time = dfsimu$T12 = dfsimu$T12.time = NA
  return(dfsimu)
}




dfsimu_creation_2times = function(n=1000, i01=0.3, i02=0.3, i12=0.3, censor_param=0.3, seed=1, ...) {
  # n=1000; i01=0.3; i02=0.3; i12=0.3; censor_param=0.3 
  dfsimu = process(n, i01, i02, i12, censor_param=NULL, seed)
  dfsimu = add_censor(dfsimu, censor_param, seed)
  head(dfsimu)
  dfsimu_creation_lastpart(dfsimu, ...)
}



dfsimu_creation_1time = function(n=1000, i01=0.3, i02=0.3, i12=0.3, censor_param=0.3, seed=1, ...) {
  # n=1000; i01=0.3; i02=0.3; i12=0.3; censor_param=0.3 
  dfsimu = process(n, i01, i02, i12, censor_param, seed)
  head(dfsimu)
  dfsimu_creation_lastpart(dfsimu, ...)
}

dfsimu_creation_lastpart = function(dfsimu, type="SIMU1", stages=c("Stage IV", "Stage II"), idx_last_indiv=0) {
  head(dfsimu)
  # dfsimu = dfsimu[,1:4]
  # colonne X1
  dfsimu = cbind(X1=1:nrow(dfsimu) + idx_last_indiv, dfsimu)
  # time in day
  dfsimu$DSS.time=dfsimu$DSS.time*365.25 
  dfsimu$PFI.time=dfsimu$PFI.time*365.25
  # adding covars
  dfsimu$type                                = type
  dfsimu$ajcc_pathologic_tumor_stage         = sample(stages                  , nrow(dfsimu), replace=TRUE)
  dfsimu$clinical_stage                      = dfsimu$ajcc_pathologic_tumor_stage
  dfsimu$gender                              = sample(c("FEMALE", "MALE")     , nrow(dfsimu), replace=TRUE)
  dfsimu$age_at_initial_pathologic_diagnosis = sample(40:90                   , nrow(dfsimu), replace=TRUE)
  return(dfsimu)
}




dfnew_creation = function(data, version="T02timePFI") {  
  # verification des colonnes numerique depuis les facteurs.
  if (!is.numeric(data$DSS.time)) {stop("DSS.time is not numeric.")}
  if (!is.numeric(data$DSS     )) {stop("DSS      is not numeric.")}
  if (!is.numeric(data$PFI.time)) {stop("PFI.time is not numeric.")}
  if (!is.numeric(data$PFI     )) {stop("PFI      is not numeric.")}
  # Change les valeurs de DSS.time et PFI.time au NA, (Car pas de sens un temps 0)
  data$DSS.time[data$DSS.time==0] <- NA #enlever les suivis nuls
  data$PFI.time[data$PFI.time==0] <- NA #enlever les suivis nuls
  # ou si PFI.time est plus grand que DSS.time car aussi pas de sens.  
  sum(data$PFI.time > data$DSS.time, na.rm=TRUE)
  # data[!is.na(data$PFI.time > data$DSS.time) & data$PFI.time > data$DSS.time,c("DSS.time", "PFI.time")]
  data[!is.na(data$PFI.time > data$DSS.time) & data$PFI.time > data$DSS.time,c("DSS.time", "PFI.time")] = NA
  # data[!is.na(data$PFI.time > data$DSS.time) & data$PFI.time > data$DSS.time,c("DSS.time", "PFI.time")]
  dim(data)
  # Creation de colonnes pour la création de dataframe de matrice de transition
  data$pfiNotDead.time <- data$PFI.time/365.25
  data$pfiNotDead <- data$PFI
  data$pfiNotDead[data$PFI.time==data$DSS.time] <- 0

  # which tiome for T02.time
  if (version=="T02timeDSS") {
    data$dssDirect.time <- data$DSS.time/365.25
  } else {
    data$dssDirect.time <- data$PFI.time/365.25    
  }

  data$dssDirect <- data$DSS
  data$dssDirect[!data$PFI.time==data$DSS.time] <- 0
  data$dssRelapse.time <- (data$DSS.time-data$PFI.time)/365.25
  data$dssRelapse.time[data$pfiNotDead==0] <- 0
  data$dssRelapse <- data$DSS
  data$dssRelapse[data$PFI.time==data$DSS.time] <- 0
  # Co variates
  # Creation de donnee simule
  data$var_nulle <-  sample(c(0,1),nrow(data), replace=TRUE)
  # data$null <-  sample(c(0,1),11160,replace=TRUE)

  # Recode les levels du stage du tumeur pour avoir moins de niveau.  
  # Fait qu'il y a seulement stage 1/2, 3/4 et NA; ou 1/2, 3/4, autre et NA  
  unique(data$ajcc_pathologic_tumor_stage)
  # "[Not Available]"        ,
  # "Stage X"                ,
  # "Stage 0"                ,
  # "[Not Applicable]"       ,
  # "[Unknown]"              ,
  # "[Discrepancy]"          ,
  stages12_labs = c(
    "I/II"   ,
    "Stage I"   ,
    "Stage IA"  ,
    "Stage IB"  ,
    "IS"        ,
    "Stage II"  ,
    "Stage IIB" ,
    "Stage IIA" ,
    "Stage IIC" ,
    "I/II NOS" 
  )
  stages34_labs = c(            
    "III/IV" ,
    "Stage III" ,
    "Stage IIIA",
    "Stage IIIB",
    "Stage IIIC",
    "Stage IV"  ,
    "Stage IVA" ,
    "Stage IVB" ,
    "Stage IVC" 
  ) 
     

  data$stage = NA
  data[data$ajcc_pathologic_tumor_stage %in% stages12_labs,]$stage = "I-II"
  data[data$ajcc_pathologic_tumor_stage %in% stages34_labs,]$stage = "III-IV"



  unique(data$clinical_stage)
  # "[Discrepancy]"
  # "[Not Available]"
  # "[Not Applicable]"
  stagebiss12_labs = c(
    "I/II"   ,
    "I"         ,      
    "Stage I"   ,      
    "Stage IA"  ,     
    "Stage IA1" ,      
    "Stage IA2" ,      
    "Stage IB"  ,     
    "Stage IB1" ,      
    "Stage IB2" ,      
    "Stage IC"  ,      
    "Stage IS"  ,      
    "IIa"       ,     
    "IIb"       ,      
    "Stage II"  ,     
    "Stage IIA" ,      
    "Stage IIA1",      
    "Stage IIA2",      
    "Stage IIB" ,      
    "Stage IIC"       
  )
  stagebiss34_labs = c(            
    "III/IV" ,
    "III"        ,     
    "Stage III"  ,     
    "Stage IIIA" ,     
    "Stage IIIB" ,     
    "Stage IIIC" ,     
    "Stage IIIC1",     
    "Stage IIIC2",     
    "IVa"        ,     
    "IVb"        ,    
    "Stage IV"   ,    
    "Stage IVA"  ,    
    "Stage IVB"  ,     
    "Stage IVC"        
  ) 
  data$stagebis = NA
  data[data$clinical_stage %in% stagebiss12_labs,]$stagebis = "I-II"
  data[data$clinical_stage %in% stagebiss34_labs,]$stagebis = "III-IV"
  table(data$stage   , useNA="ifany")
  table(data$stagebis, useNA="ifany")
  data[is.na(data$stage),]$stage = data[is.na(data$stage),]$stagebis

  data[!is.na(data$stagebis) & data$stagebis != data$stage,c("ajcc_pathologic_tumor_stage", "clinical_stage")]
  data[!is.na(data$stagebis) & data$stagebis != data$stage,"stage"] = NA
  table(data$stage   ,data$stagebis, useNA="ifany")
  

  # table(data[is.na(data$stage),]$clinical_stage)
  # table(data$clinical_stage)
  # stageBis <- as.factor(data$clinical_stage)
  # levels(stageBis) <- c("NA","NA","NA","I-II","I-II","I-II","III-IV","III-IV","III-IV",rep("I-II",14),rep("III-IV",6),"I-II",rep("III-IV",4))
  # data$stage <- sapply(1:length(stage),function(i){
  #   if(stage[i]=="NA") return(stageBis[i])
  #   else return(stage[i])
  # })
  # Creation d'un nouveau dataframe plus petit et facile a manipuler.  
  df_new <- data.frame(
    Patient=data$X1,
    DSS=data$DSS        , DSS.time=data$DSS.time/365.25 ,
    PFI=data$PFI        , PFI.time=data$PFI.time/365.25 ,
    T01=data$pfiNotDead , T01.time=data$pfiNotDead.time ,
    T02=data$dssDirect  , T02.time=data$dssDirect.time  ,
    T12=data$dssRelapse , T12.time=data$dssRelapse.time ,
    stage=data$stage                                    ,
    var_nulle=data$var_nulle                            ,
    gender=data$gender                                  ,
    age=data$age_at_initial_pathologic_diagnosis        ,
    type=data$type
  )
  # Enleve les NA
  df_new = df_new[!apply(is.na(df_new), 1, any),]
  dim(df_new)
  # Enleve les stage qui ont le niveau NA ou other
  df_new <- df_new[!is.na(df_new$stage),]
  # df_new$stage <- factor(df_new$stage,c("I-II","III-IV"))
  dim(df_new)
  df_new$age <- as.numeric(as.character(df_new$age))
  df_new = df_new[!is.na(df_new$age),] #nettoyage des NA de l'age
  # Dois les converti depuis factors pour qu'ils soit bien selectionner dans le boucle for
  # df_new <- df_new %>% mutate_at(c('type','stage','age','var_nulle','gender'), as.character)
  df_new$type       = as.character(df_new$type     )
  df_new$stage      = as.character(df_new$stage    )
  df_new$age        = as.character(df_new$age      )
  df_new$var_nulle  = as.character(df_new$var_nulle)
  df_new$gender     = as.character(df_new$gender   )
  return(df_new)
}



dfmsm_creation = function(df_new) {
  #* Section plus compliquer a comprendre
  #* Chaque personne peut atteindre max de 3 phase,  soit 1-1 (censure pour PFI et DSS), 1-2-2 (censure pour DSS), 
  #* 1-3 (pas de censure, PFI.time == DSS.time), ou 1-2-3 (pas de censure, PFI.time != DSS.time) selon leur prognosis et leur censure. 
  #* Donc il y a 3 matrice qui sont cree chaqu'un avec 4 colonne: indiv, time, state et type
  #* Dans un boucle for du longeur du dataframe des patients (qu'on n'a pas retirer) pour chaque matrix, on regarde quel phase doit placer le patient
  #* et on lui donne le bonne 'state'.  Si le personne n'a pas de troisieme phase, on met NA, qu'on enlevera apres.  
  indiv_msm <- c((1:nrow(df_new)),                                                      # tous les indiv sont dans l'état 1 (état initial)
                 (1:nrow(df_new))[df_new$T01==1],                                       # quelques indiv transit de 1 à 2
                 (1:nrow(df_new))[df_new$T12.time > 0 & df_new$T12==1],                 # quelques indiv transit de 2 à 3
                 (1:nrow(df_new))[df_new$T12.time > 0 & df_new$T12==0],                 # quelques indiv sont censurés dans l' état 2
                 (1:nrow(df_new))[df_new$T02 == 1],                                     # quelques indiv transit transit de 1 à 3
                 (1:nrow(df_new))[df_new$T02 == 0 & df_new$T01==0])                     # quelques indiv sont censurés dans l'état 1
  state_msm <- c(rep(1,nrow(df_new)),
                 rep(2,sum(df_new$T01==1)),
                 rep(3,sum(df_new$T12.time > 0 & df_new$T12==1)),
                 rep(2,sum(df_new$T12.time > 0 & df_new$T12==0)),
                 rep(3,sum(df_new$T02 == 1)),
                 rep(1,sum(df_new$T02 == 0 & df_new$T01==0)))
  time_msm <- c(rep(0,nrow(df_new)),
                df_new$T01.time[df_new$T01==1],
                (df_new$T12.time + df_new$T01.time)[df_new$T12.time > 0 & df_new$T12==1],
                (df_new$T12.time + df_new$T01.time)[df_new$T12.time > 0 & df_new$T12==0],
                df_new$T02.time[df_new$T02 == 1],
                df_new$T02.time[df_new$T02 == 0 & df_new$T01==0])
  dfmsm <- data.frame("indiv"=indiv_msm,"state"=state_msm,"time"=time_msm)
  dfmsm <- dfmsm[order(dfmsm$indiv,dfmsm$time),]
  rownames(dfmsm) <- 1:nrow(dfmsm)
  dfmsm$type <- as.character(df_new$type[dfmsm$indiv])
  # Combine les 3 matrices 
  # Ajoute les colonnes depuis df_new pour utiliser les covariants
  # Convert les donnees en bonne type de donnees
  # Enleve les donnees qui on un 'state' == NA
  dfmsm$var_nulle <- as.character(df_new[dfmsm$indiv,]$var_nulle)
  dfmsm$stage <-     as.character(df_new[dfmsm$indiv,]$stage)
  dfmsm$T01 <-       as.character(df_new[dfmsm$indiv,]$T01)
  dfmsm$T02 <-       as.character(df_new[dfmsm$indiv,]$T02)
  dfmsm$age <-       as.numeric(  df_new[dfmsm$indiv,]$age)
  # Test pour verifier les donnees sont correct
  # unique(dfmsm[dfmsm$stage=='III - IV'& dfmsm$type=='PAAD',]$indiv)
  # unique(df_new[df_new$stage=='III - IV'& df_new$type=='PAAD',]$Patient)
  # Regarde le table de transition d'etats et cree un Q matrix de base.    
  return(dfmsm)
}









get_intensity_from_data = function(d) {
  Qmat <- matrix(c(1,1,1,0,1,1,0,0,0),ncol=3,byrow = TRUE)
  dfmsm = dfmsm_creation(d) 
  Q2 <- msm::crudeinits.msm(state~time, subject=indiv, data=dfmsm, qmatrix=Qmat)
  i01 = max(0.00001, Q2[1,2])
  i02 = max(0.00001, Q2[1,3])
  i12 = max(0.00001, Q2[2,3])
  return(c(
    i01=i01, 
    i02=i02, 
    i12=i12
  ))    
}

generate_params_for_sim = function(d, kc="BRCA") {
  print(kc)
  tmp_d = d[d$type%in%kc & d$stage%in%"I-II",]
  retlo = get_intensity_from_data(tmp_d)
  nlo=nrow(tmp_d)
  names(retlo) = paste0(names(retlo), "lo")

  tmp_d = d[d$type%in%kc & d$stage%in%"III-IV",]
  rethi = get_intensity_from_data(tmp_d)
  nhi=nrow(tmp_d)
  names(rethi) = paste0(names(rethi), "hi")
  
  # tcensure < t1 / n = param obtenue a partir des données réelles (la proportion de données PFI censurée)
  # t1 censuré ~ rexp(i01+i02+i00)
  # max vraisemblance de t1 censuré ==> 1/(i01+i02+i00)
  # = mean(PFI.time)
  # i01+i02+i00 = 1/mean(PFI.time)
  # i00 = 1/mean(d$PFI.time) - (i01+i02)
  censor_param = 1/mean(d$PFI.time) - sum(get_intensity_from_data(d)[1:2])
  
  return(c(nlo=nlo, retlo, nhi=nhi, rethi, censor_param=censor_param))
}


mstate_call = function(d, cov_key="stage"){
  d$covar <- d[[cov_key]]
  tmat <- mstate::transMat(x = list(c(2, 3), c(3),c()), names = c("healthy", "ill", "dead"))
  tmat
  msebmt <- mstate::msprep(data=d, trans=tmat, time=c(NA, "T01.time", "DSS.time"), status=c(NA, "T01", "DSS"), keep="covar")
  head(msebmt)
  # msebmt[280:290,]


  # WARNING!! WARNING!! WARNING!! WARNING!! 
  msebmt$covar = d[msebmt$id,]$covar
  # WARNING!! WARNING!! WARNING!! WARNING!! 


  # msebmt[280:290,]
  # stop("probleme")

  covs <- c("covar")
  msebmt <- mstate::expand.covs(msebmt, covs, longnames = FALSE)
  head(msebmt)
  msebmt[280:290,]
  # msebmt[, c("Tstart", "Tstop", "time")] == msebmt[, c("Tstart","Tstop", "time")]
  # msebmt[, c("Tstart", "Tstop", "time")] <- msebmt[, c("Tstart","Tstop", "time")]
  library(survival)
  c_stage <- survival::coxph(survival::Surv(Tstart, Tstop, status)~strata(trans)+covar.1+covar.2+covar.3, data=msebmt, method="breslow")
  hr <- summary(c_stage)$conf.int[,-2]
  colnames(hr) <- c("HR","L","U")
  rownames(hr) <- c("mstate12", "mstate13", "mstate23")
  hr
  ret = list(hr=hr)
  return(ret)
}



run_stage_all_models = function(d) {
  
  # cox dss
  m = survival::coxph(dss~stage, d)
  sm = summary(m)
  pv_logrank = sm$sctest[3]
  tm = survival::cox.zph(m)
  tm = tm$table
  pvhz = tm[dim(tm)[1], 3]
  hr = sm$conf.int[1, 1]
  hrlb = sm$conf.int[1, 3]
  hrub = sm$conf.int[1, 4]
  res = c(hr, hrlb, hrub)
  names(res) = c("HR", "L", "U")
  res_cox_dss = res

  # cox pfi
  m = survival::coxph(pfi~stage, d)
  sm = summary(m)
  pv_logrank = sm$sctest[3]
  tm = survival::cox.zph(m)
  tm = tm$table
  pvhz = tm[dim(tm)[1], 3]
  hr = sm$conf.int[1, 1]
  hrlb = sm$conf.int[1, 3]
  hrub = sm$conf.int[1, 4]
  res = c(hr, hrlb, hrub)
  names(res) = c("HR", "L", "U")
  res_cox_pfi = res

  # msm
  Qmat = matrix(c(1,1,1,0,1,1,0,0,0), ncol=3, byrow=TRUE)
  dfmsm = dfmsm_creation(d) 
  msm_simu = mmsm(state~time, subject=dfmsm[,1], data=dfmsm, gen.inits=TRUE, qmatrix=Qmat, deathexact=3, covariates=~stage)
  hazards = msm::hazard.msm(msm_simu)[[1]]
  if(dim(hazards)[2] == 1){
    hazards = cbind(hazards, L=NA, U=NA)
  }
  rownames(hazards) = c("msm12", "msm13", "msm23")
  res_msm = hazards

  #mstate  
  res = mstate_call(d)
  res_mstate = res$hr

  ret = rbind(res_cox_dss, res_cox_pfi, res_msm, res_mstate)
  ret = data.frame(rbind(res_cox_dss, res_cox_pfi, res_msm, res_mstate))
  return(ret)
}



plot_carpet = function(d, main, ...) {
  evs = sort(unique(round(d$DSS.time,5)))
  evs = sort(unique(round(c(d$DSS.time, d$PFI.time),5)))
  nb_events = length(evs)
  ev = rep(0, nb_events)
  names(ev) = paste0("ev", evs)
  # set.seed(1)
  # d = d[sample(1:nrow(d), 300), ]
  # sum(d$T01)
  # sum(d$T02)
  # sum(d$T12)

  carpet = sapply(1:nrow(d), function(i) {
    print(i)
    ev_new = ev
    evtpfi_key = paste0("ev", round(d[i,]$PFI.time,5))
    evtdss_key = paste0("ev", round(d[i,]$DSS.time,5))

    ev_new[1:(which(names(ev_new)==evtpfi_key)-1)] = 1
    # 0, rechute ou mort ou censure
    if (d[i,]$T01) { # rechute
      if (d[i,]$T02) stop("1")
      ev_new[which(names(ev_new)==evtpfi_key):length(ev_new)] = 2
      # mort ou censure
      if (d[i,]$T12) { # mort
        ev_new[which(names(ev_new)==evtdss_key):length(ev_new)] = 3
      } else { # censure
        ev_new[which(names(ev_new)==evtdss_key):length(ev_new)] = 4
      }
    } else if (d[i,]$T02) { # mort
      if (d[i,]$T01) stop("2")  
      if (evtpfi_key!=evtdss_key) stop("3")
        ev_new[which(names(ev_new)==evtdss_key):length(ev_new)] = 3
    } else { # censure
      if (evtpfi_key!=evtdss_key) stop("4")
      ev_new[which(names(ev_new)==evtpfi_key):length(ev_new)] = 4
    }
    ev_new
  })
  colnames(carpet) = rownames(d)
  image(carpet, col=c("blue", "red", "black", "grey"), yaxt="n", xlab="events", main=main)
  axis(2, ...)
}



plot_hrs = function(results, ...) {
  par(mar=c(5.1, 10, 4.1, 2.1))
  plot(0,0,col=0, ylim=c(-nrow(results), -1), yaxt="n", xlim=c(0,10),ylab="", xlab="HR", ...)
  # min(results[,2][is.finite(results[,2])]), max(results[,3][is.finite(results[,3])])))
  abline(v=1, lty=2, col="grey")
  for (i in 1:nrow(results)) {
    col = ifelse(any(!is.finite(unlist(results[i,1:3]))), 2, 1)
    points(results[i,1], -i, pch=16, col=col)
    arrows(results[i,2],-i,results[i,3],-i, col=col, length=0)
  }
  axis(2, -(1:nrow(results)), rownames(results), las=2)
  par(mar=c(5.1, 4.1, 4.1, 2.1))  
}