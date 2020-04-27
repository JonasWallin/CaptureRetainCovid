##
# fitting maximum likelhiood usign beta-multionomial distribution
# profile likelihood to predict N
##

source('util.r')
source('stolen_function.R')
Predict.day = 10 #must be less then N
unique.p    = 5
#holidays and weekends
holidays.Sweden <- as.Date(c("2020-04-10","2020-04-11","2020-04-13"))


load("result.RData")
result_T <- result
Reported_T <- result$detected
N_T <- dim(Reported_T)[1]
Obs <- c()
for(i in 1:(N_T-Predict.day-1)){
  Obs <- c(Obs,Reported_T[i,i+Predict.day+1])
}
index_obs <- 1:N_T+Predict.day+1 <= N_T
seq_obs   <- 1:N_T
seq_obs   <- seq_obs[index_obs]
Coverage <- matrix(NA, nrow=length(seq_obs), ncol = length(seq_obs))
for(i in 5:(length(seq_obs)-1)){
  Reported <- Reported_T[1:(i+Predict.day),1:(i+Predict.day)]
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
  loglikProbBB(rep(0,2*dim(X)[2]), data_$death.remain, data_$report.new, X)$loglik
  loglik <- function(x){-loglikProbBB(x, data_$death.remain, data_$report.new, X)$loglik}
  loglik.grad <- function(x){-loglikProbBB(x, data_$death.remain, data_$report.new, X)$grad}
  res <- optim(rep(0,2*dim(X)[2]), 
               fn = loglik, 
               gr = loglik.grad,
               method="CG",
               control = list(maxit=100))
  res <- optim(res$par, 
               fn = loglik, 
               gr = loglik.grad,
               method="BFGS",
               control = list(maxit=500))
  if(res$convergence != 0){
    cat('convergence issue')
  }
  p <- dim(X)[2]
  beta_1 <- res$par[1:p]
  beta_2 <- res$par[(p+1):(2*p)]
  alpha <- exp(X%*%beta_1)
  beta <- exp(X%*%beta_2)
  Alpha_ <- matrix(NA, ncol=N, nrow=N)
  Alpha_[upper.tri(data_$report.new,diag=T)] <- alpha
  Beta_ <- matrix(NA, ncol=N, nrow=N)
  Beta_[upper.tri(data_$report.new,diag=T)] <- beta
  
  index_not <- .col(c(N,N)) - Predict.day >= .row(c(N,N))
  Alpha_[index_not] <- NA
  Beta_[index_not]  <- NA
  
  
  deaths <- death.givenParamBB(Alpha_, Beta_, Reported,Predict.day,sim=c(4000,10))
  CI <-ceiling(apply(deaths,2 , function(x){ quantile(x,c(0.025,0.975))}))
  k_ <- min(dim(deaths)[2],length(seq_obs))
  res <- CI[,(i+1):k_]-rbind(Obs[(i+1):k_],Obs[(i+1):k_])
  fig <- plot.predReport(result, CI, Predict.day)
  print(fig)
  Coverage[i,(i+1):k_] = res[1,]*res[2,]<=0 
  colnames(res) = paste(as.Date(as.numeric(colnames(res)),origin = "1970-01-01"))
  print(result$dates_report[(i+1):k_] )
  print(res)
  print(CI[2,(i+1):k_]-CI[1,(i+1):k_])
  cat('c = ',sum( res[1,]*res[2,]<=0 ),' of ',length( res[1,]*res[2,]<=0),'\n')
  a1 <- exp(beta_1[1])
  b1 <- exp(beta_2[1])  
  a1_h <- exp(beta_1[1] + beta_1[p])
  b1_h <- exp(beta_2[1]+ beta_2[p])
  cat(paste('p = ',round(a1/(a1+b1),3), ' p_h= ' ,round(a1_h/(a1_h+b1_h),3), '\n',sep=""))
  cat(paste('V[p] = ',round(100*a1*b1*(a1+b1+100)/((a1+b1)^2 * (a1+b1+1)),3),
            ' V[p_h]= ' ,round(100*a1_h*b1_h*(a1_h+b1_h+100)/((a1_h+b1_h)^2 * (a1_h+b1_h+1)),3), '\n',sep=""))
}
