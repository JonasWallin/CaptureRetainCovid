##
# fitting maximum likelhiood multinomial using logit p
# then profile likelihood of N
##
set.seed(2)
source('util.r')
source('MH.R')
source('stolen_function.R')
MCMC_sim <- 20000
burnin_p = 0.5
deaths_sim <- 10
maxusage.day = 20 #must be less then N
unique.days  = 5
true.day = 5
start.predict.day = 14# more then unique days
  
load("result.RData")
Reported = result$detected  
N <- dim(Reported)[1]

deaths_est <- apply(Reported, 1, max, na.rm=T)

data_ <- newDeaths(deaths_est,
                   Reported, 
                   maxusage.day =maxusage.day)
X <- setup_data(N, maxusage.day, result$dates_report, unique.days)

##
# seting up MCMC
##
MH_obj <- MH_setup()
MH_obj$sigma <- 0.2
MH_obj$theta <- rep(0,2*dim(X)[2])


MH_obj$Lik <- loglikProbBB

##
# mcmc loop
##
P <- matrix(NA, ncol=N, nrow=N)
Thetas <- matrix(NA, nrow=MCMC_sim, ncol = length(MH_obj$theta))
burnin <- ceiling(burnin_p*MCMC_sim)
Death_est <- matrix(NA, nrow=MCMC_sim, ncol=N)
alpha.MCMC <- rep(4, N)
p <- dim(X)[2]
Alpha <- matrix(NA, ncol=N, nrow=N)
Beta <- matrix(NA, ncol=N, nrow=N)
for(i in 1:(MCMC_sim+burnin-1)){
  if(i%%100==0){
    cat('*')
  }
  if(i%%1000==0){
    cat('t')
  }
  if(i%%10000==0){
    cat('+')
  }
  MH_obj <- MALAiter(MH_obj, TRUE,
                      death.remain = data_$death.remain,
                      report.new   = data_$report.new,
                      X            = X,
                     calcH= F)
  beta_1 <- MH_obj$theta[1:p]
  beta_2 <- MH_obj$theta[(p+1):(2*p)]
  Alpha[upper.tri(data_$report.new,diag=T)] <- exp(X%*%beta_1)
  Beta[upper.tri(data_$report.new,diag=T)]  <-  exp(X%*%beta_2)
  

  
  res <- sample.deathsBB(deaths_sim, 
                         deaths_est, 
                         Alpha, 
                         Beta, 
                         Reported, 
                         alpha.MCMC, 
                         true.day = true.day)
  deaths_est <- res$deaths
  data_ <-newDeaths(deaths_est,Reported,maxusage.day)
  
  if(i < burnin){
    alpha.MCMC[res$acc/deaths_sim > 0.3] <- alpha.MCMC[res$acc/deaths_sim > 0.3] +1
    alpha.MCMC[res$acc/deaths_sim < 0.3] <- alpha.MCMC[res$acc/deaths_sim < 0.3] -1
    alpha.MCMC[alpha.MCMC<1] <- 1
  }
  if(i >= burnin){
    Thetas[i-burnin + 1,] <-  MH_obj$theta
    Death_est[i-burnin + 1,]  <-res$deaths
  }
}
CI <-apply(Death_est,2 , function(x){ quantile(x,c(0.05,0.95))})
fig <- plot.predReport(result, CI, true.day = true.day, ymax=200)
print(fig)
ggsave(paste('data/dag_',max(result$dates),'bM.jpeg',sep=''),fig)
