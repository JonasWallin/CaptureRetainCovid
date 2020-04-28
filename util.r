require(Matrix)

##
# sample reported deaths after rep.day (if Inf) predicit day
# assuming Multionimal distribution
#
#  samples  - (int) how many samples should be done for each N (cheap)
#  deaths   - (N x 1) true number of deaths each day
#  P        - (N x N) matrix of probabilites (only upper triangular part relevant)
#  Reported - (N x N) matrix of reported deaths cumlative (only upper triangular relevant)
#  alpha    - (N x 1) stepsizes
#  true.dat  - (int) how many days after recording should one sample
##
sample.deaths <- function(samples, deaths, P, Reported, alpha, true.day = 0){

  N <- length(deaths)
  acc <- rep(0,N)
  for(i in 1:N){
    if(i > true.day){
      P_i         = P[i,i:N]
      Reported_i  = Reported[i,i:N]
      index = is.na(Reported_i)==F
      P_i         = P_i[index]
      Reported_i  = Reported_i[index]  
      lik_i <- loglikDeathsGivenProb(deaths[i], P_i, Reported_i)
      for(j in 1:samples){
        death_star <- sample((deaths[i]-alpha[i]):(deaths[i]+alpha[i]), 1)
        lik_start <- loglikDeathsGivenProb(death_star,  P_i, Reported_i)
        if(log(runif(1)) < lik_start-lik_i){
          lik_i = lik_start
          deaths[i] <- death_star
          acc[i] = acc[i] + 1
        }
      }
    }
  }
  return(list(deaths=deaths, acc = acc))
}


##
# sample reported deaths after rep.day (if Inf) predicit day using Beta binomial dist
#
#  samples  - (int) how many samples should be done for each N (cheap)
#  deaths   - (N x 1) true number of deaths each day
#  alpha     - (N x N) matrix of beta binom parameter (only upper triangular part relevant)
#  beta      - (N x N) matrix of beta binom parameter (only upper triangular part relevant)
#  Reported - (N x N) matrix of reported deaths cumlative (only upper triangular relevant)
#  alpha.MCMC    - (N x 1) stepsizes
#  true.dat  - (int) how many days after recording should one sample
##
sample.deathsBB <- function(samples, deaths, alpha, beta, Reported, alpha.MCMC, true.day = 0){
  
  N <- length(deaths)
  acc <- rep(0,N)
  for(i in 1:N){
    if(i > true.day){
      alpha_i     = alpha[i,i:N]
      beta_i      = beta[i,i:N]
      Reported_i  = Reported[i,i:N]
      index = is.na(Reported_i)==F
      alpha_i  = alpha_i[index]
      beta_i  = beta_i[index]
      Reported_i  = Reported_i[index]  
      lik_i <- loglikDeathsGivenProbBB(deaths[i],alpha_i , beta_i, Reported_i) - 0.001*deaths[i]
      for(j in 1:samples){
        death_star <- sample((deaths[i]-alpha.MCMC[i]):(deaths[i]+alpha.MCMC[i]), 1)
        lik_start <- loglikDeathsGivenProbBB(death_star, alpha_i , beta_i, Reported_i)  -0.001*death_star
        if(log(runif(1)) < lik_start-lik_i){
          lik_i = lik_start
          deaths[i] <- death_star
          acc[i] = acc[i] + 1
        }
      }
    }
  }
  return(list(deaths=deaths, acc = acc))
}

##
# fill report using binimoal beta
#
#
##
fill.ReportBB <- function(deaths, Alpha, Beta, Reported, maxusage.day){
  N <- length(deaths)
  N_2 <- dim(Alpha)[2]
  for(i in 1:N){
      Alpha_i         = Alpha[i,i:N_2]
      Beta_i          = Beta[i,i:N_2]
      Reported_i  = Reported[i,i:N_2]
      index = is.na(Reported_i)==T 
      for(j in min(which(index)):length(Reported_i)){
        p <- rbeta(1, Alpha_i[j], Beta_i[j])
        
        if(i> maxusage.day){
          Reported_i[j] <- rbinom(1, size=deaths[i] - Reported_i[j-1], prob = p) + Reported_i[j-1]
        }else{
        Reported_i[j] <- Reported_i[j-1]
        }
      }
    Reported[i,i:N_2] = Reported_i
  }
  return(Reported)
}

##
# fill report
#
#
##
fill.Report <- function(deaths, P, Reported){
  N <- length(deaths)
  N_2 <- dim(P)[2]
  for(i in 1:N){
      P_i         = P[i,i:N_2]
      Reported_i  = Reported[i,i:N_2]
      index = is.na(Reported_i)==T 
      for(j in min(which(index)):length(Reported_i)){
        Reported_i[j] <- rbinom(1, size=deaths[i] - Reported_i[j-1], prob =  P_i[j]) + Reported_i[j-1]
      }
      Reported[i,i:N_2] = Reported_i
  }
  return(Reported)
}
###
#
# posterior sampling of deaths given Prob vec
#
# P           - (N x N) probability matrix over probability of detecting reminder
# Reported    - (N x N) matrix of reported deaths cumlative (only upper triangular relevant)
# Predict.day - (int) which day to of reporting to predict
# sim         - (2 x 1) MCMC samples inner loop and outer loop
# alpha       - (N x 1) stepsizes
#
###
death.givenProb <- function(P, Reported ,Predict.day = Inf,sim=c(2000,10), alpha = NULL){

  N <- dim(P)[1]
  if(is.null(alpha))
    alpha <- rep(4, N)

  deaths <- matrix(NA, nrow=sim[1], ncol = N)
  deaths_est <- apply(Reported, 1, max, na.rm=T)

  burnin = ceiling(0.3*sim[1])
  for(i in 1:(sim[1] + burnin -1)){

    res <- sample.deaths(sim[2],deaths_est, P, Reported, alpha,rep.day=Predict.day)
    deaths_est <- res$deaths
    if(i < burnin){
      alpha[res$acc/sim[2] > 0.3] <- alpha[res$acc/sim[2] > 0.3] +1
      alpha[res$acc/sim[2] < 0.3] <- alpha[res$acc/sim[2] < 0.3] -1
      alpha[alpha<1] <- 1
    }else{
      deaths[i - burnin +1,] = deaths_est
    }
  }
  return(deaths)
}




##
# log liklihood of obseving report given death and prob
#  deaths  - (int) true number of deaths
#  p       - (n x 1) probability of report
#  report  - (n x 1) reported deaths culmative each date
##
loglikDeathsGivenProb <- function(death, p, report){

  if(death < max(report,na.rm=T))
    return(-Inf)
  n <- length(p)
  if(is.na(report[1])){
    #we dont have data from day one

    ndeaths <- diff(report[is.na(report)==F])
    report_adj <- report[1:(n-1)]
    report_adj <- report_adj[is.na(report_adj)==F]
  }else{
    ndeaths <- c(report[1],diff(report))
    if(n>1){
      report_adj <- c(0, report[1:(n-1)])
    }else{
      report_adj <- 0
    }
  }

  ndeaths[ndeaths <0 ] = 0

  return(sum(dbinom(ndeaths,  death - report_adj,  prob = p[is.na(p)==F], log=T )))
}

##
# build holiday covariates vector for holiday
#
##
#  holidays - (N x 1) true if day is holiday false else
##
buildXholiday <- function(N,holidays){
  ##
  # base matrix
  ##
  index_base <- t(matrix(rep(holidays,N),ncol = N))

  index_base <- index_base[upper.tri(index_base,diag=T)]
  index_base <- 1* index_base
  #sparse matrix
  i_base <- 1:length(index_base)
  i_base <- i_base[index_base==1]
  return(sparseMatrix(j=rep(1,length(i_base)),i=i_base,  dims=c(length(index_base), 1)))
}

buildXall <- function(N){
  ##
  # base matrix
  ##
  index_base <- t(matrix(1:N^2,ncol = N, nrow=N))

  index_base <- index_base[upper.tri(index_base,diag=T)]
  #sparse matrix
  i_base <- 1:length(index_base)
  j_base <- 1:length(index_base)
  return(sparseMatrix(j=j_base,i=j_base))
}

##
# build day effect matrix
#  nDayEffects - number of speical days effect (1- first day, 2- first + second day)
#  N - number of days
#  nDayEffects - how many day covariate effect to create
##
buildXdayeffect <- function(N, nDayEffects = 1){
  nDayEffects <- min(N,nDayEffects)
  index_days <- toeplitz(c( (1:nDayEffects), rep(0,N -nDayEffects)))
  index_days <- index_days[upper.tri(index_days,diag=T)]
  j_ <-  index_days[index_days>0]
  i_base <- 1:length(index_days)
  i_ <-  i_base[index_days>0]
  return(sparseMatrix(i=i_,j=j_ , dims=c(length(index_days), nDayEffects)  ) )
}
##
# build day mixed effect matrix
#  nDayEffects - number of speical days effect (1- first day, 2- first + second day)
#  N - number of days
#  nDayEffects - how many day covariate effect to create
##
buildXmixeddayeffect <- function(N, nDayEffects = 1){

  index = rep(0,N)
  index[nDayEffects] = 1
  index_days <- toeplitz(index)
  index_days <- index_days[upper.tri(index_days,diag=T)]
  j_ <-  1:sum(index_days[index_days>0])
  i_base <- 1:length(index_days)
  i_ <-  i_base[index_days>0]
  return(sparseMatrix(i=i_,j=j_ , dims=c(length(index_days), sum(index_days[index_days>0])))  )
}
##
# building an X matrix such that each day has a fixed effect
#
##
buildXday <- function(N){
  ##
  # base matrix
  ##
  index_base <- t(matrix(rep(1:N,N),ncol = N))
  index_base <- index_base[upper.tri(index_base,diag=T)]
  #sparse matrix
  j_base <- 1:length(index_base)
  #adding mean effect
  j_ <- c(index_base, rep(N+1,length(index_base)))
  i_ <- c(j_base, 1:length(index_base))
  #day effects
  if(nDayEffects > 0){
    index_days <- toeplitz(c( N+1 + (1:nDayEffects), rep(0,N -nDayEffects)))
    index_days <- index_days[upper.tri(index_days,diag=T)]
    j_ <- c(j_, index_days[index_days>0])
    i_ <- c(i_, j_base[index_days>0])
  }
  return(sparseMatrix(i=i_,j=j_))
}
##
# transforms data,
# we remove data if negative..
#
#  deaths  - (N x 1) how many has died (thruth)
#  reports - (N x N) reported cumlative deaths
#  maxusage.day - (int) only use data up to maxusage.days i.e.
#                       reports[i,i + maxusage.day - 1]
#
##
newDeaths <-function(deaths, reports,maxusage.day = -1){
  newreport <- reports
  diff.report <- t(diff(t(reports)))
  newreport[upper.tri(newreport)] <- diff.report[upper.tri(diff.report,T)]
  newreport[newreport<0 & is.na(newreport)==F]=0 #fake
  death.rem <- diag(deaths)
  dr<-deaths-reports
  N <- length(deaths)
  dr <- dr[1:(N-1),1:(N-1)]
  death.rem[upper.tri(newreport)] <- dr[upper.tri(dr,T)]
  death.rem[lower.tri(death.rem)] <- NA
  diag(death.rem)[is.na(diag(reports))] <- NA
  if(maxusage.day>0){
    for(i in 1:length(death.rem)){
      if(i + maxusage.day - 1 < N ){
        death.rem[i,(i+maxusage.day):N] = NA
        newreport[i,(i+maxusage.day):N] = NA
      }
    }
  }
  #removing small bugs in reporting
  newreport[death.rem <0 & is.na(death.rem)==F] = 0
  death.rem[death.rem <0& is.na(death.rem)==F] = 0
  index <- (death.rem< newreport) &  is.na(death.rem) ==F
  newreport[index] =death.rem[index] 
  return(list(death.remain = death.rem, report.new = newreport))
}


##
# likelihood observations given deaths and probabilites
# model is:
# newreport[upper.tri] \sim   Bin(deaths[upper.tri], logit(X\beta))
#
#' @return obj - list
#'               $loglik  - logliklihood
#'               $grad    - gradient
#'               $hessian - Hessian of the likelihood
##
loglikProb <- function(beta, death.remain, report.new, X){

  N <- dim(report.new)[1]
  
  index <- upper.tri(death.remain,diag = T)
  y = report.new[index]
  n = death.remain[index]
  index = is.na(y)==F
  X <- X[index,]
  p = as.vector(1/(1+exp(-X%*%beta)))
  y <- y[index]
  n <- n[index]
  lik <- sum(dbinom(y,size = n, prob = as.vector(p), log = T))

  grad <- as.vector(t(X)%*%as.vector(y-n*p))

  Hessian <- -t(X)%*%diag(as.vector(n*p*(1-p)))%*%X
  Hessian <- diag(diag(Hessian))

  return(list(loglik = lik, grad = grad, Hessian = Hessian))
}
loglikPrior <- function(beta, n, sigma, res){
  n <- length(sigma)
  res$grad[1:n] <-  res$grad[1:n]-beta[1:n]/(sigma^2)
  diag(res$Hessian) <- diag(res$Hessian) - c(1/sigma^2, rep(0,length(beta)-n))
  res$loglik <- res$loglik -sum(beta[1:n]^2/(2*sigma^2))
  return(res)
}

loglik <- function(beta,death.remain, report.new, X, sigma){
  res <- loglikProb(beta, death.remain, report.new, X)
  res <- loglikPrior(beta, n, sigma, res)
  return(res)
}


#####
##
# beta binomial part
##
#####



##
# likelihood observations given deaths and probabilites
# model is:
# newreport[upper.tri] \sim   BB(deaths[upper.tri], exp(X\beta_1), exp(X\beta_2))
#
#' @return obj - list
#'               $loglik  - logliklihood
#'               $grad    - gradient
##
loglikProbBB <- function(beta, death.remain, report.new, X, calcH=T){
  
  p <- dim(X)[2]
  beta_1 <- beta[1:p]
  beta_2 <- beta[(p+1):(2*p)]
  N <- dim(report.new)[1]
  
  index <- upper.tri(death.remain,diag = T)
  y = report.new[index]
  n = death.remain[index]
  index = is.na(y)==F
  X <- X[index,]
  alpha <- exp(X%*%beta_1)
  beta  <- exp(X%*%beta_2) #unfortunate name
  if(min(1/(alpha+beta)) <0.001)
    return(list(loglik = -Inf, grad = 0, Hessian = 0))
  y <- y[index]
  n <- n[index]
  
  lik <- sum(dBB(y, size = n, alpha = alpha, beta = beta, log.p=T))
  if(is.na(lik))
    return(list(loglik = -Inf, grad = 0, Hessian = 0))
  
  grad_lik <- grad_dBB(y, size = n, alpha = alpha, beta = beta)
  grad_alpha <- as.vector(t(X)%*%(alpha*grad_lik$grad_alpha))
  grad_beta  <- as.vector(t(X)%*%(beta*grad_lik$grad_beta))
  if(calcH){
    H_ <- Hessian_dBB(y, size = n, alpha = alpha, beta = beta)
    H_12 <-  t(X)%*%diag(as.vector(alpha*beta*H_$grad_alpha_beta))%*%X
    H_11 <-  t(X)%*%diag(as.vector(alpha^2*H_$grad_alpha_alpha))%*%X
    H_22 <-  t(X)%*%diag(as.vector(beta^2*H_$grad_beta_beta))%*%X
    Hessian <- rbind(cbind(H_11, H_12) ,
                     cbind(H_12, H_22) )
    Hessian  <- diag(diag(Hessian) )
    diag(Hessian) <- diag(Hessian)
  }else{
    Hessian =NULL
  }
  grad =  c(grad_alpha,grad_beta)
    if(sum(is.na(grad))>0)
      return(list(loglik = -Inf, grad = 0, Hessian = 0))
  return(list(loglik = lik, grad = grad, Hessian = Hessian))
}
##
# PRobability of beta binomial
##
dBB<- function(x, size, alpha, beta, log.p = F){
  const = lgamma(size + 1) - lgamma(x + 1) -
    lgamma(size - x + 1)
  ld <- const + 
    lgamma(x + alpha) + lgamma(size - x + beta) -
    lgamma(size + alpha + beta) + 
    lgamma(alpha + beta) - 
    lgamma(alpha) - lgamma(beta)
  if(log.p==F){
    return(exp(ld))
  }
  return(ld)
}
##
# gradient of log of probability beta binomial
# return list, 
# [[1]] grad alpha
# [[2]] grad beta
##
grad_dBB<- function(x, size, alpha, beta){
  grad_alpha <- digamma(x + alpha) -
    digamma(size + alpha + beta) + 
    digamma(alpha + beta) -
    digamma(alpha)
  grad_beta <- digamma(size - x + beta) -
    digamma(size + alpha + beta) + 
    digamma(alpha + beta) -
    digamma(beta)
  return(list(grad_alpha = grad_alpha,
              grad_beta  = grad_beta ))
}
##
# Hessian of log of probability beta binomial
# return list, 
# [[1]] grad alpha
# [[2]] grad beta
##
Hessian_dBB<- function(x, size, alpha, beta){
  grad_alpha_alpha <- trigamma(x + alpha) -
    trigamma(size + alpha + beta) + 
    trigamma(alpha + beta) -
    trigamma(alpha)
  grad_beta_beta <- trigamma(size - x + beta) -
    trigamma(size + alpha + beta) + 
    trigamma(alpha + beta) -
    trigamma(beta)
  
  grad_alpha_beta <- -trigamma(size + alpha + beta) + 
    trigamma(alpha + beta) 
  return(list(grad_alpha_alpha = grad_alpha_alpha,
              grad_beta_beta  = grad_beta_beta,
              grad_alpha_beta = grad_alpha_beta))
}

##
# log liklihood of obseving report given death and prob
# density is Beta binomial
#  deaths  - (int) true number of deaths
#  alpha       - (n x 1) bb parameter 1
#  beta        - (n x 1) bb parameter 2
#  report  - (n x 1) reported deaths culmative each date
##
loglikDeathsGivenProbBB <- function(death, alpha, beta, report){
  
  if(death < max(report,na.rm=T))
    return(-Inf)
  n <- length(alpha)
  if(is.na(report[1])){
    #we dont have data from day one
    
    ndeaths <- diff(report[is.na(report)==F])
    report_adj <- report[1:(n-1)]
    report_adj <- report_adj[is.na(report_adj)==F]
  }else{
    ndeaths <- c(report[1],diff(report))
    if(n>1){
      report_adj <- c(0, report[1:(n-1)])
    }else{
      report_adj <- 0
    }
  }
  
  ndeaths[ndeaths <0 ] = 0
  
  return(sum(dBB(x =ndeaths, 
                 size = death - report_adj, 
                 alpha = alpha[is.na(alpha)==F], 
                 beta = beta[is.na(alpha)==F],
                 log.p=T)))
}


###
#
# posterior sampling of deaths given Prob vec
#
#  alpha     - (N x N) matrix of beta binom parameter (only upper triangular part relevant)
#  beta      - (N x N) matrix of beta binom parameter (only upper triangular part relevant)
# Reported    - (N x N) matrix of reported deaths cumlative (only upper triangular relevant)
# Predict.day - (int) which day to of reporting to predict
# sim         - (2 x 1) MCMC samples inner loop and outer loop
# alpha       - (N x 1) stepsizes
#
###
death.givenParamBB <- function(alpha, beta, Reported,Predict.day = Inf,sim=c(2000,10), alpha.MCMC = NULL){
  
  N <- dim(alpha)[1]
  if(is.null(alpha.MCMC))
    alpha.MCMC <- rep(4, N)
  
  deaths <- matrix(NA, nrow=sim[1], ncol = N)
  deaths_est <- apply(Reported, 1, max, na.rm=T)
  
  burnin = ceiling(0.3*sim[1])
  for(i in 1:(sim[1] + burnin -1)){
    res <- sample.deathsBB(sim[2],deaths_est, alpha, beta, Reported, alpha.MCMC,rep.day=Predict.day)
    deaths_est <- res$deaths
    if(i < burnin){
      alpha.MCMC[res$acc/sim[2] > 0.3] <- alpha.MCMC[res$acc/sim[2] > 0.3] +1
      alpha.MCMC[res$acc/sim[2] < 0.3] <- alpha.MCMC[res$acc/sim[2] < 0.3] -1
      alpha.MCMC[alpha.MCMC<1] <- 1
    }else{
      deaths[i - burnin +1,] = deaths_est
    }
  }
  return(deaths)
}

##
# buidling the X matrix
#
#
# unique.prob - (int) how many of days shall we have unique probabilites
##
setup_data <- function(N, Predict.day, dates_report, unique.prob=NULL){
  holidays.Sweden <- as.Date(c("2020-04-10","2020-04-13"))
  X <- buildXdayeffect(N,Predict.day)
  #holidays and weekends
  
  result$dates_report <- as.Date(dates_report)
  holidays <- weekdays(dates_report)%in%c("Sunday","Saturday") | c(dates_report)%in%c(holidays.Sweden)
  
  holidays.tommorow <- weekdays(dates_report + 1)%in%c("Sunday","Saturday") |
    (result$dates_report+1)%in%c(holidays.Sweden)
  Xhol <- buildXholiday(N,holidays)
  #Xhol.tommrow <- buildXholiday(N,holidays.tommorow)
  if(is.null(unique.prob)==F)
    X <- cbind(X[,1:unique.prob],rowSums(X[,1:unique.prob])==0)
  
  X <- cbind(X, Xhol)
  return(X)
}


