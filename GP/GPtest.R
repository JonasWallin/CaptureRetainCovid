###
# testing the GP code
###
source('../MH.R')
source('GPutil.R')
library(rstan)
set.seed(9)
sim <- 50000
n <- 30
tau  = 10
theta <- log(100)+cumsum(rnorm(n))/tau
lambda <- exp(theta)#apply(as.matrix(theta),1,function(x){max(x,0.1)})

N <- rpois(n, lambda)


MH_obj <- MH_setup()
MH_obj$sigma <- 0.1
MH_obj$theta <- log(N)



MH_obj$Lik <- function(x, N, tau) {
  prior_list <- prior_s1d(x, tau)
  lik_list   <- lik_poisson(N, x)
  return(list(loglik = lik_list$loglik + prior_list$loglik, 
              grad = lik_list$grad  + prior_list$grad, 
              Hessian = lik_list$Hessian + prior_list$Hessian,
              tau = prior_list$tau))
}


thetaM <- matrix(0, nrow=sim, ncol = n)
n <- length(theta)
L <- toeplitz(c(-1,1, rep(0,n-1)))
L[lower.tri(L,diag=F)] <- 0
L <- L[-(n+1),]
L <- L[,-1]
L<-as(L, "sparseMatrix")
A = solve(sqrt(tau)*L)
tau_vec <- rep(0,sim)
for(i in 1:sim){
  MH_obj <- MALAiter(MH_obj, TRUE,
                     N = N,
                     tau = tau) 
  thetaM[i,] <-  as.vector(MH_obj$theta)
  tau_vec[i] <- MH_obj$res$tau
  tau        <- MH_obj$res$tau
}
thetaM <- thetaM[ceiling(0.5*sim):sim,]
par(mfrow=c(1,2))
plot(thetaM[,1])
plot(colMeans(thetaM), ylim=c(min(colMeans(thetaM)),1.1*max(theta)))
lines(theta,col='red')

data.list <- list(N=N,t = length(N),
                  tau = tau) 
stan_out = stan(file = "carStan.stan",
                data= data.list,
                iter = 5000,
                warmup = ceiling(0.5*5000),
                chains = 2)
Theta_stan <- extract(stan_out)$theta
points(colMeans(Theta_stan),col='blue')