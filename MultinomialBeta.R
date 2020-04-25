##
# fitting maximum likelhiood usign beta-multionomial distribution
# profile likelihood to predict N
##

source('util.r')
source('stolen_function.R')
Predict.day = 14 #must be less then N


load("result.RData")
Reported <- result$detected
N <- dim(Reported)[1]
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
#holidays and weekends
holidays.Sweden <- as.Date(c("2020-04-10","2020-04-13"))
result$dates_report <- as.Date(result$dates_report)
holidays <- weekdays(result$dates_report)%in%c("Sunday","Saturday") |c(result$dates_report)%in%c(holidays.Sweden)

holidays.tommorow <- weekdays(result$dates_report + 1)%in%c("Sunday","Saturday") |
                (result$dates_report+1)%in%c(holidays.Sweden)
Xhol <- buildXholiday(holidays)
Xhol.tommrow <- buildXholiday(holidays.tommorow)
X <- cbind(X, Xhol,  Xhol.tommrow[,1]*X[,dim(X)[2]])
loglikProbBB(rep(0,2*dim(X)[2]), data_$death.remain, data_$report.new, X)$loglik
loglik <- function(x){-loglikProbBB(x, data_$death.remain, data_$report.new, X)$loglik}
loglik.grad <- function(x){-loglikProbBB(x, data_$death.remain, data_$report.new, X)$grad}
res <- optim(rep(0,2*dim(X)[2]), 
             fn = loglik, 
             gr = loglik.grad,
             method="BFGS",
             control = list(maxit=500))
if(res$convergence != 0){
  cat('convergence issue')
}
#res$value
# 268.6968
# res$value holidays still matters
# 241.4545

#not true prob but mean
p <- dim(X)[2]
beta_1 <- res$par[1:p]
beta_2 <- res$par[(p+1):(2*p)]
alpha <- exp(X%*%beta_1)
beta <- exp(X%*%beta_2)
P <- matrix(NA, ncol=N, nrow=N)
Pvec <- alpha/(alpha+beta)
P[upper.tri(data_$report.new,diag=T)] <- Pvec

Alpha <- matrix(NA, ncol=N, nrow=N)
Alpha[upper.tri(data_$report.new,diag=T)] <- alpha
Beta <- matrix(NA, ncol=N, nrow=N)
Beta[upper.tri(data_$report.new,diag=T)] <- beta

index_not <- .col(c(N,N)) - Predict.day >= .row(c(N,N))
Alpha[index_not] <- NA
Beta[index_not]  <- NA
P[index_not]     <- NA


deaths <- death.givenParamBB(Alpha, Beta, Reported,Predict.day,sim=c(2000,10))
CI <-apply(deaths,2 , function(x){ quantile(x,c(0.05,0.95))})

fig <- plot.predReport(result, CI, Predict.day)
print(fig)

print(fig)
ggsave(paste('data/dag_',Predict.day+1,'_',max(result$dates),'bM.jpeg',sep=''),fig)
