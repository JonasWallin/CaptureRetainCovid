##
# fitting maximum likelhiood multinomial using logit p
# then profile likelihood of N
##

source('util.r')
source('stolen_function.R')
Predict.day = 13 #must be less then N


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
holidays <- weekdays(result$dates_report)%in%c("Sunday","Saturday") |result$dates_report%in%c(holidays.Sweden)

holidays.tommorow <- weekdays(result$dates_report + 1)%in%c("Sunday","Saturday") |
  (result$dates_report+1)%in%c(holidays.Sweden)
Xhol <- buildXholiday(holidays)
Xhol.tommrow <- buildXholiday(holidays.tommorow)
X <- cbind(X, Xhol,  Xhol.tommrow[,1]*X[,dim(X)[2]])
loglikProb(rep(0,dim(X)[2]), data_$death.remain, data_$report.new, X)
loglik <- function(x){-loglikProb(x, data_$death.remain, data_$report.new, X)$loglik}
loglik.grad <- function(x){-loglikProb(x, data_$death.remain, data_$report.new, X)$grad}
res <- optim(rep(0,dim(X)[2]), fn = loglik, gr = loglik.grad, method="CG")
P <- matrix(NA, ncol=N, nrow=N)
#res$value
#   385.5662
#res$value holidays matter
# 308.2527
Pvec <- as.vector(1/(1+ exp(-X%*%res$par)))
P[upper.tri(data_$report.new,diag=T)] <- Pvec

deaths <- death.givenProb(P = P,
                          Reported = result$detected,
                          Predict.day = Predict.day,
                          sim=c(1000,10))


CI <-apply(deaths,2 , function(x){ quantile(x,c(0.05,0.95))})
fig <- plot.predReport(result, CI, Predict.day)
print(fig)
ggsave(paste('data/dag_',Predict.day+1,'_',max(result$dates),'Ml.jpeg',sep=''),fig)