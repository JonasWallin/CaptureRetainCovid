##
# fitting maximum likelhiood usign beta-multionomial distribution
# profile likelihood to predict N
##

source('util.r')
source('stolen_function.R')
Predict.day = 10 #must be less then N
unique.p    = 5
#holidays and weekends
holidays.Sweden <- as.Date(c("2020-04-10","2020-04-13"))


load("result.RData")
result_T <- result
Reported_T <- result$detected
N_T <- dim(Reported_T)[1]
Obs <- c()
for(i in 1:(N_T-Predict.day)){
  Obs <- c(Obs,Reported_T[i,i+Predict.day])
}
index_obs <- 1:N_T+Predict.day <= N_T
seq_obs   <- 1:N_T
seq_obs   <- seq_obs[index_obs]
Coverage <- matrix(NA, nrow=length(seq_obs), ncol = length(seq_obs))
for(i in 3:(length(seq_obs)-1)){
  Reported <- Reported_T[1:(i+Predict.day-1),1:(i+Predict.day-1)]
  N <- dim(Reported)[1]
  result <- result_T
  result$detected <- Reported
  result$dates <- result$dates[1:N]
  result$dates_report <- result$dates_report[1:N]
  data_report <- result$dates_report[1:N]
  deaths_est <- apply(Reported, 1, max, na.rm=T)

  data_ <- newDeaths(deaths_est,
                     Reported, 
                     rep.day =Predict.day)
  #remove the "future" data
  data_$death.remain[(N - Predict.day+1):N,] = NA
  data_$report.new[(N - Predict.day+1):N,] = NA
  ##
  # building covariate matrix
  ##
  X <- buildXdayeffect(N,Predict.day)
  result$dates_report <- as.Date(data_report)
  holidays <- weekdays(data_report)%in%c("Sunday","Saturday") |c(data_report)%in%c(holidays.Sweden)
  
  holidays.tommorow <- weekdays(data_report + 1)%in%c("Sunday","Saturday") |
    (data_report+1)%in%c(holidays.Sweden)
  Xhol <- buildXholiday(holidays)
  Xhol.tommrow <- buildXholiday(holidays.tommorow)
  X <- cbind(X[,1:unique.p],rowSums(X[,(unique.p+1):dim(X)[2]]))
  X <- cbind(X, Xhol )
  X<- X[,colSums(X)>2]
  loglikProb(rep(0,dim(X)[2]), data_$death.remain, data_$report.new, X)
  loglik <- function(x){-loglikProb(x, data_$death.remain, data_$report.new, X)$loglik}
  loglik.grad <- function(x){-loglikProb(x, data_$death.remain, data_$report.new, X)$grad}
  res <- optim(rep(0,dim(X)[2]), fn = loglik, gr = loglik.grad, method="CG")
  P <- matrix(NA, ncol=N, nrow=N)
  Pvec <- as.vector(1/(1+ exp(-X%*%res$par)))
  P[upper.tri(data_$report.new,diag=T)] <- Pvec
  
  deaths <- death.givenProb(P = P,
                            Reported = result$detected,
                            Predict.day = Predict.day,
                            sim=c(1000,10))
  CI <-ceiling(apply(deaths,2 , function(x){ quantile(x,c(0.05,0.95))}))
  k_ <- min(dim(deaths)[2],length(seq_obs))
  res_CI <- CI[,(i+1):k_]-rbind(Obs[(i+1):k_],Obs[(i+1):k_])
  fig <- plot.predReport(result, CI, Predict.day)
  print(fig)
  Coverage[i,(i+1):k_] = res_CI[1,]*res_CI[2,]<=0 
  print(res_CI)
  cat('c = ',sum( res_CI[1,]*res_CI[2,]<=0 ),' of ',length( res_CI[1,]*res_CI[2,]<=0),'\n')
  print(round(1/(1+exp(-res$par[1:(unique.p+1)])),2))
  print(round(1/(1+exp(-res$par[1] - res$par[(unique.p+2)])),2))
}
