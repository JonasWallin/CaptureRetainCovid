##
# fitting maximum likelhiood usign beta-multionomial distribution
# profile likelihood to predict N
##

source('util.r')
source('stolen_function.R')
Predict.day = 10 #predicit up to day 
unique.p    = 5


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
X <- setup_data(N, Predict.day, data_$dates_report)
##
# optim
##
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
