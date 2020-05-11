##
# fitting maximum likelhiood multinomial using logit p
# then profile likelihood of N
##
set.seed(2)
source('util.r')
source('MH.R')
source('stolen_function.R')
source('./GP/GPutil.R')
MCMC_sim <- 40000
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
MH_obj$sigma <- 0.1
MH_obj$theta <- rep(0,2*dim(X)[2])


MH_obj$Lik <- loglikProbBB

##
# seting up MCMC GP
##
MH_obj_GP <- MH_setup()
MH_obj_GP$sigma <- 0.1
#probably better do something smarter
MH_obj_GP$theta <- log(deaths_est+1) 
tau <- 0.1
MH_obj_GP$Lik <- function(x, N, tau) {
  prior_list <- prior_s1d(x, tau)
  lik_list   <- lik_poisson(N, x)
  return(list(loglik = lik_list$loglik + prior_list$loglik, 
              grad = lik_list$grad  + prior_list$grad, 
              Hessian = lik_list$Hessian + prior_list$Hessian,
              tau = prior_list$tau))
}

##
# mcmc loop
##
P <- matrix(NA, ncol=N, nrow=N)
Thetas <- matrix(NA, nrow=MCMC_sim, ncol = length(MH_obj$theta))

theta_GP <-  matrix(NA, nrow=MCMC_sim, ncol = length(MH_obj_GP$theta))
tau_vec <- rep(0,MCMC_sim)
burnin <- ceiling(burnin_p*MCMC_sim)
Death_est <- matrix(NA, nrow=MCMC_sim, ncol=N)
alpha.MCMC <- rep(4, N)
p <- dim(X)[2]
Alpha <- matrix(NA, ncol=N, nrow=N)
Beta <- matrix(NA, ncol=N, nrow=N)
X_next <- setup_data(N+1, 
                     maxusage.day, 
                     c(result$dates_report,result$dates_report[length(result$dates_report)]+1), 
                     unique.days)
Reported_fill <- cbind(Reported, matrix(NA, nrow=N,ncol = 1))
Alpha_next <- matrix(NA, N + 1,N + 1)
Beta_next  <- matrix(NA, N + 1,N + 1)
pred_next <- matrix(NA, nrow=MCMC_sim, ncol=N)
for(i in 1:(MCMC_sim+burnin-1)){
  if(i%%100==0){
    cat('*')
  }
  if(i%%1000==0){
    cat('t')
  }
  if(i%%10000==0){
    cat(' ',i/10000,' ')
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
  

  prior <- function(N, i){ N *MH_obj_GP$theta[i] - lgamma(N+1)}
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
  
  
  MH_obj_GP <- MALAiter(MH_obj_GP, TRUE,
                     N = deaths_est,
                     tau = tau) 
  tau        <- MH_obj_GP$res$tau
  
  if(i < burnin){
    alpha.MCMC[res$acc/deaths_sim > 0.3] <- alpha.MCMC[res$acc/deaths_sim > 0.3] +1
    alpha.MCMC[res$acc/deaths_sim < 0.3] <- alpha.MCMC[res$acc/deaths_sim < 0.3] -1
    alpha.MCMC[alpha.MCMC<1] <- 1
  }
  if(i >= burnin){
    Thetas[i-burnin + 1,] <-  MH_obj$theta
    Death_est[i-burnin + 1,]  <-res$deaths
    
    theta_GP[i-burnin + 1,] <-  as.vector(MH_obj_GP$theta)
    tau_vec[i-burnin + 1] <- MH_obj_GP$res$tau
    Alpha_next[upper.tri(Alpha_next,diag=T)] <- exp(X_next%*%beta_1)
    Beta_next[upper.tri(Beta_next,diag=T)]   <- exp(X_next%*%beta_2)
    Reported_sample <-fill.ReportBB(res$deaths,  
                                    Alpha_next[1:dim(Reported_fill)[1],1:dim(Reported_fill)[2]],
                                    Beta_next[1:dim(Reported_fill)[1],1:dim(Reported_fill)[2]], 
                                    Reported_fill,
                                    maxusage.day = maxusage.day)
    pred_next[i-burnin + 1,] <-Reported_sample[,dim(Reported)[2]+1]
    
  }
}
CI <-apply(Death_est,2 , function(x){ quantile(x,c(0.05,0.95))})
fig <- plot.predReport(result, CI, true.day = true.day, ymax=min(max(CI)+5,200))
print(fig)


##
# addd seven day rolling average
##
MA <- rep(0, N)
MaxR <-  apply(Reported,1, max, na.rm=T)
for(i in 1:N){
  MA_temp <- 
  MA[i] <- mean(MaxR[max(1,i-6):i])
}
roll_average = data.frame( date =result$dates[1:(N-7)], 
                           Reported = MA[1:(N-7)] ) 
roll_average$cumReported <- cumsum(roll_average$Reported)
fig <- plot.predReport(result, CI, true.day = true.day, ymax=min(max(CI)+5,200))
fig <- fig  + geom_line(data= roll_average,
                        mapping=aes(y = Reported, x = date),
                        inherit.aes = FALSE,
                        lwd=1.,
                        color='blue')
fig <- fig  + geom_line(data= data.frame(date =result$dates, 
                                         int = colMeans(exp(theta_GP))),
                        mapping=aes(y = int, x = date),
                        inherit.aes = FALSE,
                        lwd=1.,
                        color='black')
print(fig)

ggsave(paste('data/dag_',max(result$dates),'bMGP.jpeg',sep=''),fig)
jpeg(paste('data/dag_',max(result$dates),'_hist1','bMGP.jpeg',sep=''))
par(mfrow=c(1,2))
hist(Death_est[,dim(Death_est)[2]-1],
     xlim=c(0,200),
     breaks=100, 
     probability = T, 
     main=paste('döda ',result$dates[dim(Death_est)[2]-1],sep=""))
hist(Death_est[,dim(Death_est)[2]-2],
     xlim=c(0,200),
     breaks=100, 
     probability = T,
     main=paste('döda ',result$dates[dim(Death_est)[2]-2],sep=""))
dev.off()

CI_pred <-apply(pred_next,2 , function(x){ quantile(x,c(0.05,0.95))})
CI_pred[1,1:true.day]=deaths_est[1:true.day]
CI_pred[2,1:true.day]=deaths_est[1:true.day]

fig <- plot.predReport(result, CI_pred, true.day = true.day, ymax=200)
fig <- fig +  ggtitle(paste("prediktion för rapport ",result$dates[length(result$dates)]+1,sep=""))
ggsave(paste('data/rapport_',max(result$dates)+1,'bMGP.jpeg',sep=''),fig)
print(data.frame(date= result$dates, CI_l = t(CI_pred)[,1], CI_u = t(CI_pred)[,2]))

Q <-apply(apply(Death_est,1,cumsum),1,quantile, probs=c(0.05,0.5,0.95))

cumpred <- data.frame(dates = result$dates, deaths = Q[2,], lQ = Q[1,], uQ= Q[3,])
default_theme = set_default_theme()
fig1 <-ggplot(cumpred) + 
  default_theme+
  geom_line(aes(y = deaths, x = dates)) +
  geom_ribbon(aes(x=dates,ymin = lQ, ymax = uQ), alpha = .2) +
  scale_x_date(date_breaks = "4 day",
               expand = c(0, 0),
               limits = as.Date(c('2020-04-20',max(cumpred$dates))))

print(fig1)
ggsave(paste('data/cumlativedeaths_',max(result$dates),'bMGP.jpeg',sep=''),fig1)

MAsim <- matrix(NA, nrow= MCMC_sim, ncol = N)
for(i in 1:N){
  MAsim[ ,i] <- apply(Death_est,1,function(x){mean(x[max(1,i-6):i])})
  MAsim[ ,i] <- apply(Death_est,1,function(x){mean(x[max(1,i-6):i])})
}
QMA <-apply(MAsim,2,quantile, probs=c(0.05,0.5,0.95))
MApred <- data.frame(date = result$dates, deaths = QMA[2,], lQ = QMA[1,], uQ= QMA[3,])
fig <- ggplot(data= roll_average,
              mapping=aes(y = Reported, x = date)) +geom_line()  + 
            scale_x_date(date_breaks = "4 day",
                         expand = c(0, 0),
                         limits = as.Date(c(min(cumpred$dates),max(cumpred$dates)-7))) +
          ylim(c(70,110)) +
          geom_line(data= MApred,
                      mapping=aes(y = deaths, x = date),
                      color='blue') + 
          geom_line(data= MApred,
                    mapping=aes(y = lQ, x = date),
                    linetype = "dashed",
                    color='blue') + 
          geom_line(data= MApred,
                    mapping=aes(y = uQ, x = date),
                    linetype = "dashed",
                    color='blue') 
  print(fig)
  ggsave(paste('data/movingaverage_',max(result$dates),'bMGP.jpeg',sep=''),fig)
  
###
# mix check
###   
x11()
par(mfrow=c(2,1))
plot(exp(theta_GP[,N]))
plot(tau_vec)
save(list = ls(all.names = TRUE), file = "GP_prior.RData", envir = .GlobalEnv)