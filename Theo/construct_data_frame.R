construct_data_frame <- function(df,choice.type,choice.end,covariate=NULL){
  endpoints <- c("OS","DSS","PFI","DFI")
  endpoints.time <- c("OS.time","DSS.time","PFI.time","DFI.time")
  type <- c()
  time <- c()
  survVal <- c()
  survValUp <- c()
  survValDo <- c()
  covar <- c()
  end <- c()
  for(j in choice.end){
    for(i in choice.type){
      if(!is.null(covariate)){
        for(l in unique(covariate)){
          if(!all(is.na(df[[endpoints.time[[j]]]][df$type==i & covariate==l])) & !all(is.na(df[[endpoints[[j]]]][df$type==i & covariate==l]))){
            y <- Surv(df[[endpoints.time[[j]]]][df$type==i & covariate==l],as.numeric(df[[endpoints[[j]]]][df$type==i & covariate==l]))
            sf <- survfit(y~1)
            survVal <- c(survVal,sf$surv)
            survValUp <- c(survValUp,sf$upper)
            survValDo <- c(survValDo,sf$lower)
            time <- c(time,sf$time)
            type <- c(type,rep(i,length(sf$time)))
            covar <- c(covar,rep(l,length(sf$time)))
            end <- c(end,rep(endpoints[[j]],length(sf$time)))
          }
        }
      }
      else{
          if(!all(is.na(df[[endpoints.time[[j]]]][df$type==i & covariate==l])) & !all(is.na(df[[endpoints[[j]]]][df$type==i & covariate==l]))){
            y <- Surv(df[[endpoints.time[[j]]]][df$type==i & covariate==l],as.numeric(df[[endpoints[[j]]]][df$type==i & covariate==l]))
            sf <- survfit(y~1)
            survVal <- c(survVal,sf$surv)
            survValUp <- c(survValUp,sf$upper)
            survValDo <- c(survValDo,sf$lower)
            time <- c(time,sf$time)
            type <- c(type,rep(i,length(sf$time)))
            covar <- c(covar,rep(l,length(sf$time)))
            end <- c(end,rep(endpoints[[j]],length(sf$time)))
          }
      }
    }
  }
  time <- time/365.25
  dfSurv <- data.frame(time,survVal,survValUp,survValDo,type,end,covar)
  return(dfSurv)
}
