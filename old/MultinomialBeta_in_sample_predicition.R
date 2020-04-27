##
# fitting maximum likelhiood usign beta-multionomial distribution
# profile likelihood to predict N
##
source('MultinomialBeta.R')
Obs <- c()
for(i in 1:(N-Predict.day)){
  Obs <- c(Obs,Reported[i,i+Predict.day])
}
index_obs <- 1:N+Predict.day <=N
Coverage <- c()
for(i in 2:Predict.day){
  cat('i = ',i,' of ',Predict.day,'\n')
  index <- .col(c(N,N)) - i + 1>= .row(c(N,N))
  
  Reported_fake <- Reported
  Reported_fake[index] = NA
  deaths <- death.givenParamBB(Alpha, Beta, Reported_fake,Predict.day,sim=c(4000,10))
  CI <-apply(deaths,2 , function(x){ quantile(x,c(0.05,0.95))})
  res <- CI[,index_obs]-rbind(Obs,Obs)
  Coverage <- rbind(Coverage,(res[1,]*res[2,]<=0))
  fig <- plot.predReport(result, CI, Predict.day)
  print(fig)
  print(CI[2,]-CI[1,])
  cat('c = ',sum( res[1,]*res[2,]<=0 ),' of ',length( res[1,]*res[2,]<=0),'\n')
}
