##
# fitting maximum likelhiood multinomial using logit p
# then profile likelihood of N
##
source('util.r')
source('MH.R')
source('stolen_function.R')
MCMC_sim <- 4000
burnin_p = 0.2
deaths_sim <- 10
maxusage.day = 20 #must be less then N
unique.days  = 5
true.day = 5
start.predict.day = 30# more then unique days
  
load("result.RData")
Reported_T = result$detected  
N_T <- dim(Reported_T)[1]

deaths_est_T <- apply(Reported_T, 1, max, na.rm=T)

data_T <- newDeaths(deaths_est_T,
                   Reported_T, 
                   maxusage.day =maxusage.day)
X_T <- setup_data(N_T, maxusage.day, result$dates_report, unique.days)
predicition.list <-list()
for(j in start.predict.day:(N_T-1)){
  cat('\nj = ',j,'\n')
  load("result.RData")
  Reported <- result$detected[1:j,1:j]
  dates_report <- result$dates_report[1:j]
  N <- dim(Reported)[1]
  deaths_est[1:true.day] = deaths_est_T[1:true.day] #days with known deaths
  deaths_est <- apply(Reported, 1, max, na.rm=T)
  
  data_ <- newDeaths(deaths_est,
                     Reported, 
                     maxusage.day =maxusage.day)
  
  
  ##
  # building covariate matrix
  ##
  X <- setup_data(N, maxusage.day,  dates_report, unique.days)
  
  ##
  # seting up MCMC
  ##
  MH_obj <- MH_setup()
  MH_obj$sigma <- 0.1
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
  pred_set <-array(NA, dim = c(j,N_T-j,MCMC_sim)) 
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
                           true.day = true.day,
                           prior=prior)
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
      
      Reported_fill <- cbind(Reported, matrix(NA, nrow=j,ncol = (N_T-j)))
      Alpha_T <- matrix(NA, N_T,N_T)
      Beta_T  <- matrix(NA, N_T,N_T)
      Alpha_T[upper.tri(Alpha_T,diag=T)] <- exp(X_T%*%beta_1)
      Beta_T[upper.tri(Beta_T,diag=T)]   <- exp(X_T%*%beta_2)
      Alpha_T <- Alpha_T[1:dim(Reported_fill)[1],1:dim(Reported_fill)[2]]
      Beta_T  <- Beta_T[1:dim(Reported_fill)[1],1:dim(Reported_fill)[2]]
      Reported_sample <-fill.ReportBB(res$deaths, 
                                      Alpha_T,
                                      Beta_T, 
                                      Reported_fill, 
                                      maxusage.day = maxusage.day)
      Reported_sample[,dim(Reported)[2]+1]
      pred_set[,,i-burnin + 1] <- Reported_sample[,(j+1):dim(Reported_sample)[2]]
    }
  }
  
  predicition.list[[j]]<- list(samples = pred_set, 
                               truth = Reported_T[1:j,(j+1):N_T],
                               j = j,
                               thetas = Thetas,
                               MH_obj = MH_obj,
                               deaths = Death_est)
  #traceplots
  par(mfrow=c(1,2))
  plot(Thetas[,1])
  plot(Death_est[,N])
  
  
  CI <-apply(Death_est,2 , function(x){ quantile(x,c(0.05,0.95))})
  result_j <- result
  result_j$detected <- result$detected[1:j,1:j]
  result_j$dates <- result$dates[1:j]
  result_j$dates <- result$dates_report[1:j]
  fig <- plot.predReport(result_j, CI)
  print(fig)

}
coverage <- matrix(NA, nrow = N_T-1, ncol = N_T-start.predict.day)
for(j in 1:(N_T-2)){
  if(is.null(predicition.list[[j]])==F){
    Quan<-apply(predicition.list[[j]]$samples,c(1,2), quantile,probs=c(0.025,0.975))
    for(k in 1:dim(Quan)[3]){
      truth <- as.matrix(predicition.list[[j]]$truth)[,k]
      res <- Quan[,,k] - rbind(truth,truth)
      res <- res[,-(1:true.day)]
      coverage[j,k] <- mean(res[1,]*res[2,]<=0)
    }
  }
}