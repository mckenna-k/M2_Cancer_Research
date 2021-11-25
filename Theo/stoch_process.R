S <- 1:3
# process <- function(S,t,mat){
#   if(t >1){
#     print(t)
#     return(process(sample(S,1,prob=mat[S,]),t-1,mat))
#   }
#   else return(S)
# }
# process(S,Time,mat)
n <- 10
S <- 1:3
s <- 1
Time <- 1:1000

mat <- function(S, t, t0) {
  return(matrix(
    c(
      atan(1000 - (t - t0)) / pi + 0.5,
      (0.5 - atan(1000 - (t - t0)) / pi) * 0.9,
      (0.5 - atan(1000 - (t - t0)) / pi) * 0.1,
      0,
      atan(1000 - (t - t0)) / pi + 0.5,
      0.5 - atan(1000 - (t - t0)) / pi,
      0,
      0,
      1
    ),
    byrow = TRUE,
    ncol = length(S),
    nrow = length(S)
  ))
}

cont <- c(0,0,0)
plot(0:1000,rep(0,1001),type="l",ylim=c(0,3))
for(j in 1:n){
  res <- c(s)
  compt <- 1
  t0 <- 0
  while(res[[length(res)]] != S[[length(S)]]){
    if(length(res) < 1) if(res[[length(res)]] > res[[length(res)-1]]) t0 <- compt
    res <- c(res,sample(S,1,prob=mat(S,compt,t0)[res[[length(res)]],]))
    compt <- compt + 1
  }
  points(cumsum(as.vector(table(res))),1:length(unique(res)),type="b",col=j)
  # res <- c(res,c(1,2,3))
  # cont <- cont + table(res)
}
# cont <- cont - n
# sum(cont)
# cont
