###
# capture retian regular Maxmimum likelihood
# then profile likelihood of N
###

source('util.r')
source('stolen_function.R')
Predict.day = 13 #must be less then N


load("result.RData")
N <- dim(result$detected)[1]
p <- rep(0,Predict.day-1)
N_sum <- rep(0,Predict.day)
for(i in 1:Predict.day)
{
  for(j in 1:(N - Predict.day+1)){
    N_sum[i] <-  N_sum[i]  + result$detected[j ,j + i - 1]
  }
}
yobs <- diff(c(0,N_sum[1:(Predict.day-1)]))
Nobs <- (N_sum[Predict.day]-c(0,N_sum[1:(Predict.day-2)]))
p <- yobs/Nobs
llik <- sum(yobs *log(p) + (Nobs-yobs) * log(1-p))
P <- toeplitz(c(p,rep(NA,N-length(p))))
P[lower.tri(P)] <- NA

Reported <- result$detected
Pdeaths <- matrix(NA, nrow=N, ncol=200)

deaths <- death.givenProb(P = P,
                          Reported = result$detected,
                          Predict.day = Predict.day,
                          sim=c(1000,10))


CI <-apply(deaths,2 , function(x){ quantile(x,c(0.05,0.95))})
fig <- plot.predReport(result, CI, Predict.day)
print(fig)

print(fig)
ggsave(paste('data/dag_',Predict.day+1,'_',max(result$dates),'Msimple.jpeg',sep=''),fig)