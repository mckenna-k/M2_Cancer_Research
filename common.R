process <- function(n, i01, i02, i12, censor1=NULL, censor2=NULL, seed=1){
  set.seed(seed)
  if(is.null(censor1) | is.null(censor2)){
    t1 <- rexp(n,i01+i02) #temps de passage de 0 à 1 ou 0 à 2
    s1 <- rbinom(n,1,i02/(i01+i02))+1 #détermination de l'état atteint (1 ou 2)
    t2 <- rexp(n,i12) #temps de passage de l'état 1 à 2 à partir du temps d'arrivée à 1
    df <- data.frame("DSS"=rep(1,n),"DSS.time"=ifelse(s1==1,t1+t2,t1),
                     "PFI"=rep(1,n),"PFI.time"=t1,
                     "T01"=-(s1-2),"T01.time"=t1,
                     "T02"=s1-1,"T02.time"=ifelse(s1==2,t1,t1+t2),
                     "T12"=-(s1-2),"T12.time"=ifelse(s1==1,t2,0))
  }  else {
    t1 <- rexp(n,i01+i02+censor1) #temps de passage de 0 à 1 ou 0 à 2
    s1 <- apply(rmultinom(n, size = 1, prob = c(censor1/(i01+i02+censor1),i01/(i01+i02+censor1),i02/(i01+i02+censor1)))==1,2,which)-1 #détermination de l'état atteint (1 ou 2)
    t2 <- rexp(n,i12+censor2)
    s2 <- rbinom(n,1,i12/(i12+censor2))
    df <- data.frame("DSS"=ifelse(s1 > 0,1,0),"DSS.time"=ifelse(s1==1,t1+t2,t1),
                     "PFI"=ifelse(s1 > 0,1,0),"PFI.time"=t1,
                     "T01"=ifelse(s1==1,1,0),"T01.time"=t1,
                     "T02"=ifelse(s1-1 > 0,1,0),"T02.time"=ifelse(s1==2 | s1==0,t1,t1+t2),
                     "T12"=ifelse(s1==1 & s2==1,1,0),"T12.time"=ifelse(s1==1,t2,0))
  }
  return(df)
}



generate_params_for_sim = function(kc="BRCA", df_new) {
  Qmat <- matrix(c(1,1,1,0,1,1,0,0,0),ncol=3,byrow = TRUE)
  df_msm = df_new_2_df_msm(df_new)
  Q2 <- crudeinits.msm(state~time, subject=indiv,data = df_msm[df_msm$type == kc,],qmatrix= Qmat)
  i01 = Q2[1,2]
  i02 = Q2[1,3]
  i12 = Q2[2,3]
  # d = df_msm[df_msm$type == kc,]
  # head(d)
  # table(d$T01)
  # censor1 = exp(-)
  # censor2 = exp(-)
  
  return(c(i01=i01, i02=i02, i12=i12))  

}







df_new_2_df_msm = function(df_new) {
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
  df_msm <- data.frame("indiv"=indiv_msm,"state"=state_msm,"time"=time_msm)
  df_msm <- df_msm[order(df_msm$indiv,df_msm$time),]
  rownames(df_msm) <- 1:nrow(df_msm)
  df_msm$type <- as.character(df_new$type[df_msm$indiv])
  # Combine les 3 matrices 
  # Ajoute les colonnes depuis df_new pour utiliser les covariants
  # Convert les donnees en bonne type de donnees
  # Enleve les donnees qui on un 'state' == NA
  df_msm$var_nulle <- as.character(df_new[df_msm$indiv,]$var_nulle)
  df_msm$stage <-     as.character(df_new[df_msm$indiv,]$stage)
  df_msm$T01 <-       as.character(df_new[df_msm$indiv,]$T01)
  df_msm$T02 <-       as.character(df_new[df_msm$indiv,]$T02)
  df_msm$age <-       as.numeric(  df_new[df_msm$indiv,]$age)
  # Test pour verifier les donnees sont correct
  # unique(df_msm[df_msm$stage=='III - IV'& df_msm$type=='PAAD',]$indiv)
  # unique(df_new[df_new$stage=='III - IV'& df_new$type=='PAAD',]$Patient)
  # Regarde le table de transition d'etats et cree un Q matrix de base.    
  return(df_msm)
}


