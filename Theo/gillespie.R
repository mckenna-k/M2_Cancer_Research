library(survival)
library(ggplot2)

#génération d'une trajectoire
process <- function(mat,s,end){
  time <- c()
  state <- c()
  new_state <- s
  new_time <- 0
  while(new_time < end){
    time <- c(time,new_time)
    state <- c(state,new_state)
    lambda <- sum(mat[new_state,])
    if(lambda == 0) break
    new_time <- rexp(1,lambda) + new_time
    if(new_time > end) break
    trans_not_null <- mat[new_state,] > 0
    vec <- c(0,cumsum(mat[new_state,][trans_not_null]))
    var_unif <- runif(1,0,lambda)
    i <- 1
    while(new_state == state[[length(state)]]){
      if(var_unif < vec[[i]] | vec[[i+1]] <= var_unif) i <- i + 1
      else if(var_unif == vec[[1]]) new_state <- which(trans_not_null)[[1]]
      else new_state <- which(trans_not_null)[[i]]
    }
  }
  return(data.frame(time,state))
}

#graphe de plusieurs trajectoires
plot.process <- function(mat,s,end,B){
  df <- process(mat,s,end)
  num_curve <- rep(1,nrow(df))
  for(i in 2:B){
    new_df <- process(mat,s,end)
    df <- rbind(df,new_df,c(end,new_df$state[[nrow(new_df)]]))
    num_curve <- c(num_curve,rep(i,nrow(new_df)+1))
  }
  df <- cbind(df,"num_curve"=as.factor(num_curve))
  ggplot(df)+
    geom_step(aes(time,state,color=num_curve))+
    theme(panel.background=element_blank(),panel.grid = element_line(size=0.1,color = "grey"),legend.position = "None")+
    scale_y_continuous(breaks = sort(unique(df$state),decreasing = FALSE))
}

#courbe de survie et de risque cumulé
plot.process.KM <- function(mat,s,end,B){
  time <- c()
  event <- c()
  for(i in 1:B){
    df <- process(mat,s,end)
    if(2 %in% df$state){
      time <- c(time,df$time[df$state==2])
      event <- c(event,1)
      if(3 %in% df$state){
        time <- c(time,df$time[df$state==3],end)
        event <- c(event,1,0)
      }
      else{
        time <- c(time,end,end)
        event <- c(event,0,0)
      }
    }
    else{
      time <- c(time,end,end)
      event <- c(event,0,0)
      if(3 %in% df$state){
        time <- c(time,df$time[df$state==3])
        event <- c(event,1)
      }
      else{
        time <- c(time,end)
        event <- c(event,0)
      }
    }
  }
  transi <- rep(c("1 -> 2","2 -> 3","1 -> 3"),B)
  df <- data.frame(time,event,transi)
  support <- survfit(Surv(df$time,df$event)~transi,data=df)
  new_df <- data.frame("time"=rep(support$time,2),"val"=c(support$surv,support$cumhaz),"transitions"=rep(rep(names(support$strata),support$strata),2),"type"=rep(c("Survival","Cumulative hazard"),each=length(support$time)))
  
  ggplot(new_df)+
    geom_step(aes(time,val,color=transitions))+
    facet_wrap(~type,scales = "free_y")+
    theme_minimal()
}

mat <- matrix(c(0,0.3,0.1,0.2,0,0.1,0,0,0),ncol=3,nrow = 3,byrow=TRUE)
mat <- matrix(c(0,0.3,0.1,0,0,0.1,0,0,0),ncol=3,nrow = 3,byrow=TRUE)
end <- 20
s <- 1 
B <- 30
plot.process(mat,s,end,B)
plot.process.KM(mat,s,end,B)