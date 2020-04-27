##
# fitting maximum likelhiood multinomial using logit p
# then profile likelihood of N
##

source('util.r')
source('MH.R')
source('stolen_function.R')
MCMC_sim <- 1000
deaths_sim <- 10
maxusage.day = 13 #must be less then N
unique.days  = 5
true.day = 5 #treat N of the first five days as the final thruth

load("result.RData")
Reported <- result$detected
N <- dim(Reported)[1]
deaths_est <- apply(Reported, 1, max, na.rm=T)

data_ <- newDeaths(deaths_est,
                   Reported, 
                   maxusage.day =maxusage.day)


##
# building covariate matrix
##
X <- setup_data(N, maxusage.day, result$dates_report, unique.days)

##
# seting up MCMC
##
MH_obj <- MH_setup()
MH_obj$sigma <- 0.5
MH_obj$theta <- rep(0,dim(X)[2])


MH_obj$Lik <- loglikProb
theta_mean <- rep(0,dim(X)[2])
thetas <- c()

##
# mcmc loop
##
P <- matrix(NA, ncol=N, nrow=N)
Thetas <- matrix(NA, nrow=MCMC_sim, ncol = length(MH_obj$theta))
burnin <- ceiling(0.1*MCMC_sim)
Death_est <- matrix(NA, nrow=MCMC_sim, ncol=N)
alpha <- rep(4, N)

for(i in 1:(MCMC_sim+burnin-1)){
  
  MH_obj <- MALAHiter(MH_obj, TRUE,
                      death.remain = data_$death.remain,
                      report.new   = data_$report.new,
                      X            = X)

    P[upper.tri(data_$report.new,diag=T)] <- as.vector(1/(1+ exp(-X%*%MH_obj$theta)))
    res <- sample.deaths(deaths_sim ,deaths_est, P, Reported, alpha, true.day = true.day)
    
    deaths_est <- res$deaths
    data_ <-newDeaths(deaths_est,Reported,maxusage.day)
    
    if(i < burnin){
      alpha[res$acc/deaths_sim > 0.3] <- alpha[res$acc/deaths_sim > 0.3] +1
      alpha[res$acc/deaths_sim < 0.3] <- alpha[res$acc/deaths_sim < 0.3] -1
      alpha[alpha<1] <- 1
    }
    if(i >= burnin){
      Thetas[i-burnin + 1,] <-  MH_obj$theta
      Death_est[i-burnin + 1,]  <-res$deaths
    }
}

#traceplots
par(mfrow=c(1,2))
plot(Thetas[,1])
plot(Death_est[,20])


CI <-apply(Death_est,2 , function(x){ quantile(x,c(0.05,0.95))})
fig <- plot.predReport(result, CI)
print(fig)
