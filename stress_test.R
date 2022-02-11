process <- function(n, i01, i02, i12, i00=NULL, i11=NULL, seed=1){
  set.seed(seed)
  if(is.null(i00) | is.null(i11)){
    t1 <- rexp(n,i01+i02) #temps de passage de 0 à 1 ou 0 à 2
    s1 <- rbinom(n,1,i02/(i01+i02))+1 #détermination de l'état atteint (1 ou 2)
    t2 <- rexp(n,i12) #temps de passage de l'état 1 à 2 à partir du temps d'arrivée à 1
    df <- data.frame("DSS"=rep(1,n),"DSS.time"=ifelse(s1==1,t1+t2,t1),
                     "PFI"=rep(1,n),"PFI.time"=t1,
                     "T01"=-(s1-2),"T01.time"=t1,
                     "T02"=s1-1,"T02.time"=ifelse(s1==2,t1,t1+t2),
                     "T12"=-(s1-2),"T12.time"=ifelse(s1==1,t2,0))
  }  else {
    t1 <- rexp(n,i01+i02+i00) #temps de passage de 0 à 1 ou 0 à 2 ou censor 
    s1 <- apply(rmultinom(n, size = 1, prob = c(i00/(i01+i02+i00), i01/(i01+i02+i00), i02/(i01+i02+i00)))==1,2,which)-1 #détermination de l'état atteint (1 ou 2 ou 0)
    t2 <- rexp(n,i12+i11) # idem
    s2 <- rbinom(n,1,i12/(i12+i11)) # 2 ou 1
    df <- data.frame("DSS"=ifelse(s1 > 0,1,0),"DSS.time"=ifelse(s1==1,t1+t2,t1),
                     "PFI"=ifelse(s1 > 0,1,0),"PFI.time"=t1,
                     "T01"=ifelse(s1==1,1,0),"T01.time"=t1,
                     "T02"=ifelse(s1-1 > 0,1,0),"T02.time"=ifelse(s1==2 | s1==0,t1,t1+t2),
                     "T12"=ifelse(s1==1 & s2==1,1,0),"T12.time"=ifelse(s1==1,t2,0))
  }
  return(df)
}




#Import les données
data = openxlsx::read.xlsx("TCGA-CDR-SupplementalTableS1.xlsx", sheet = "TCGA-CDR", na.strings="#N/A")
head(data)
data = data[,c("X1", "DSS", "DSS.time", "PFI", "PFI.time", "type", "ajcc_pathologic_tumor_stage", "clinical_stage", "gender",  "age_at_initial_pathologic_diagnosis")]
head(data)




df_simu1 <- process(n=1000,0.3,0.3,0.3,0.3,0.3,97)
df_simu1_toadd = df_simu1[,1:4]
df_simu1_toadd = cbind(X1=1:nrow(df_simu1_toadd) + nrow(data), df_simu1_toadd)
head(df_simu1_toadd)

df_simu1_toadd$DSS.time=df_simu1_toadd$DSS.time*365.25 
df_simu1_toadd$PFI.time=df_simu1_toadd$PFI.time*365.25

df_simu1_toadd$type                                = "SIMU1"
df_simu1_toadd$ajcc_pathologic_tumor_stage         = sample(c("Stage IV", "Stage II")        , nrow(df_simu1_toadd), replace=TRUE)
df_simu1_toadd$clinical_stage                      = df_simu1_toadd$ajcc_pathologic_tumor_stage
df_simu1_toadd$gender                              = sample(c("FEMALE", "MALE")                    , nrow(df_simu1_toadd), replace=TRUE)
df_simu1_toadd$age_at_initial_pathologic_diagnosis = sample(na.omit(data$age_at_initial_pathologic_diagnosis), nrow(df_simu1_toadd))

colnames(data)
colnames(df_simu1_toadd)

data = rbind(data, df_simu1_toadd)


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
data[!is.na(data$PFI.time > data$DSS.time) & data$PFI.time > data$DSS.time,c("DSS.time", "PFI.time")]
data[!is.na(data$PFI.time > data$DSS.time) & data$PFI.time > data$DSS.time,c("DSS.time", "PFI.time")] = NA
data[!is.na(data$PFI.time > data$DSS.time) & data$PFI.time > data$DSS.time,c("DSS.time", "PFI.time")]
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
stage <- as.factor(data$ajcc_pathologic_tumor_stage)
levels(stage) <- c("NA","NA","NA","NA","I-II","I-II","other","I-II","I-II","I-II","I-II","I-II","I-II","I-II","III-IV","III-IV","III-IV","III-IV","III-IV","III-IV","III-IV","III-IV","other")
stageBis <- as.factor(data$clinical_stage)
levels(stageBis) <- c("NA","NA","NA","I-II","I-II","I-II","III-IV","III-IV","III-IV",rep("I-II",14),rep("III-IV",6),"I-II",rep("III-IV",4))
data$stage <- sapply(1:length(stage),function(i){
  if(stage[i]=="NA") return(stageBis[i])
  else return(stage[i])
})





# Creation d'un nouveau dataframe plus petit et facile a manipuler.  
df_new <- data.frame(
  Patient=data$X,
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
df_new <- df_new[-which(is.element(df_new$stage,c("NA","other"))),]
df_new$stage <- factor(df_new$stage,c("I-II","III-IV"))
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




head(df_new)

head(df_simu1)

head(df_new[df_new$type=="SIMU1",colnames(df_simu1)])

sum(df_new[df_new$type=="SIMU1",]$T01       != df_simu1$T01      )
sum(df_new[df_new$type=="SIMU1",]$T02       != df_simu1$T02)      
sum(df_new[df_new$type=="SIMU1",]$T12       != df_simu1$T12 )     
sort(df_new[df_new$type=="SIMU1",]$T01.time  - df_simu1$T01.time)
sort(df_new[df_new$type=="SIMU1",]$T02.time  - df_simu1$T02.time )
sort(df_new[df_new$type=="SIMU1",]$T12.time  - df_simu1$T12.time )


