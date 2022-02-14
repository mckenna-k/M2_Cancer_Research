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


process <- function(n, i01, i02, i12, censor_param=NULL, seed=1){
  set.seed(seed)
  t1 <- rexp(n,i01+i02) #temps de passage de 0 à 1 ou 0 à 2
  s1 <- rbinom(n,1,i02/(i01+i02))+1 #détermination de l'état atteint (1 ou 2)
  t2 <- rexp(n,i12) #temps de passage de l'état 1 à 2 à partir du temps d'arrivée à 1
  if(is.null(censor_param)){
    dfsimu <- data.frame("DSS"=rep(1,n),"DSS.time"=ifelse(s1==1,t1+t2,t1),
                     "PFI"=rep(1,n),"PFI.time"=t1,
                     "T01"=-(s1-2),"T01.time"=t1,
                     "T02"=s1-1,"T02.time"=ifelse(s1==2,t1,t1+t2),
                     "T12"=-(s1-2),"T12.time"=ifelse(s1==1,t2,0))
  }  else {
    set.seed(seed*2+1)
    tcensure <- rexp(n,censor_param) #temps de censure pour chaque individu
    print(paste0("nb_censored_DSS", sum(t1+t2 > tcensure)))
    print(paste0("nb_censored_PFI", sum(t1 > tcensure)))
    tmin0102 = t1
    t1 <- apply(matrix(c(tcensure,tmin0102),byrow = FALSE,ncol=2),1,min) #min temps de passage de 0 à 1 et de 0 à 2 et du temps de censure
    s1 <- ifelse(tmin0102 < tcensure,s1,0) #détermination de l'état atteint (1 ou 2 ou censure)
    tmin12 <- t2
    t2 <- apply(matrix(c(tmin12,tcensure-t1),byrow = FALSE,ncol=2),1,min) #min temps de passage de 1 à 2 et du temps de censure
    s2 <- ifelse(tmin12 < tcensure,1,0) #détermination de la censure de la transition 1 -> 2
    dfsimu <- data.frame("DSS"=ifelse(s1 == 2 | (s1 == 1 & s2 == 1),1,0),"DSS.time"=ifelse(s1==1,t1+t2,t1),
                     "PFI"=ifelse(s1 > 0,1,0),"PFI.time"=t1,
                     "T01"=ifelse(s1==1,1,0),"T01.time"=t1,
                     "T02"=ifelse(s1-1 > 0,1,0),"T02.time"=ifelse(s1==2 | s1==0,t1,t1+t2),
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
  
  print(paste0("nb_censored_DSS", sum(dfsimu$DSS.time > tcensure)))
  print(paste0("nb_censored_PFI", sum(dfsimu$PFI.time > tcensure)))

  dfsimu$DSS      = ifelse(dfsimu$DSS.time > tcensure, 0       , 1)
  dfsimu$DSS.time = ifelse(dfsimu$DSS.time > tcensure, tcensure, dfsimu$DSS.time)
  dfsimu$PFI      = ifelse(dfsimu$PFI.time > tcensure, 0       , 1)
  dfsimu$PFI.time = ifelse(dfsimu$PFI.time > tcensure, tcensure, dfsimu$PFI.time)

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




dfnew_creation = function(data) {  
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
  data$dssDirect.time <- data$DSS.time/365.25
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









generate_params_for_sim = function(kc="BRCA", df_new) {
  # kc = "BRCA"
  Qmat <- matrix(c(1,1,1,0,1,1,0,0,0),ncol=3,byrow = TRUE)
  dfmsm = dfmsm_creation(df_new)
  dfmsm = dfmsm[dfmsm$type == kc,]
  n = nrow(dfmsm)  
  Q2 <- msm::crudeinits.msm(state~time, subject=indiv, data=dfmsm[dfmsm$type == kc,],qmatrix= Qmat)
  i01 = Q2[1,2]
  i02 = Q2[1,3]
  i12 = Q2[2,3]
  d = df_new[df_new$type == kc,]
  propcensor_param = 1 - sum(d$T01 + d$T02) / nrow(d)
  propi11 = 1 - sum(d$T12) / sum(d$T01)
  censor_param = propcensor_param * (i01 + i02) / (1 - propcensor_param)
  i11 = propi11*i12 / (1-propi11)
  # 1 / mean(d$PFI) - i01 -i02 # likelihodd based (Theo)
   
  # head(d)
  # table(d$T01)
  # table(d$T02)
  # table(d$T12)
  # censor_param = exp(-1)
  # i11 = exp(-1)
  return(c(n=n, i01=i01, i02=i02, i12=i12, censor_param=censor_param, i11=i11))  
}







