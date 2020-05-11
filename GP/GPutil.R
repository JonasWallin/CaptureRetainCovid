library(Matrix)

#'
#' second order random walk model (intrisinct)
#'  d theta /dt = N(0,\tau)
#'  @param theta - (nx1)
#'  @param tau   - precision
#'  @param alpha - prior parameter
#'  @param beta  - prior parameter
prior_s1d <-function(theta, tau, alpha = 0.01, beta = 0.01){
  
  n <- length(theta)
  L <- toeplitz(c(-1,1, rep(0,n-2)))
  L[lower.tri(L,diag=F)] <- 0
  L <- L[-n,]
  L<-as(L, "sparseMatrix")
  lik <- -tau/2 * sum((L%*%theta)^2)
  Q = t(L)%*%L
  grad <- -tau * (Q%*%theta)
  Hessian <- - tau *Q
  ##
  # sample tau
  ##
  sigma   <- 1/rgamma(1, shape= (n-1)/2+1 + alpha, rate = -lik/tau + beta)
  tau_sample <- 1/sigma
  return(list(loglik = lik, grad = grad, Hessian = Hessian, tau=tau_sample))
}



#'
#' Poisson log likelihood
#' N     - (n x 1) Poisson observation
#' theta - (n x 1) log intens
#' 
lik_poisson<- function(N, theta){
  lambda <- exp(theta)
  lik <- sum(N * theta -  lambda)
  grad <- N  - lambda
  Hessian <- - SparseM::diag(lambda)
  return(list(loglik = lik, grad = grad, Hessian = Hessian))
}
#'
#' Poisson log likelihood
#' N     - (n x 1) Poisson observation
#' theta - (m x 1) log intens
#' A     - (n x m) link matrix
#' 
lik_poisson_noise<- function(N, theta,A){
  llambda <- A%*%theta
  lambda <- exp(llambda)
  lik <- sum(N * (llambda) -  lambda)
  grad <- t(A)%*%(N  - lambda)
  Hessian <- - t(A)%*%SparseM::diag(lambda)%*%A
  return(list(loglik = lik, grad = grad, Hessian = Hessian))
}
#'
#' Poisson log likelihood
#' N     - (n x 1) Poisson observation
#' theta - (m x 1) log intens
#' A     - (n x m) link matrix
#' 
lik_poisson_noise<- function(N, theta,A){
  llambda <- as.vector(A%*%theta)
  lambda <- as.vector(exp(llambda))
  lik <- sum(N * (llambda) -  lambda)
  grad <- t(A)%*%(N  - lambda)
  Hessian <- - t(A)%*%SparseM::diag(lambda)%*%A
  return(list(loglik = lik, grad = grad, Hessian = Hessian))
}

#'
#' second order random walk model (intrisinct)
#'  d theta /dt = N(0,\tau)
#'  @param theta - (nx1)
#'  @param tau   - precision
prior_noise <-function(theta){
  
  n <- length(theta)
  lik <- - 0.5 * sum(theta[2:n]^2)
  grad <- - c(0,theta[2:n])
  Hessian <- - SparseM::diag(n)
  Hessian[1,1] <- 0
  return(list(loglik = lik, grad = grad, Hessian = Hessian))
}



#'
#' Poisson log likelihood
#' N     - (n x 1) Poisson observation
#' theta - (n x 1) log intens
#' 
lik_poisson2<- function(N, theta){
  lambda <- apply(as.matrix(theta),1,function(x){max(x,0.1)})
  lik <- sum(N * log(lambda) -  lambda)
  grad <- (theta>0.1) * (N/lambda - 1)
  Hessian <- - SparseM::diag((theta>0.1)/lambda^2)
  return(list(loglik = lik, grad = grad, Hessian = Hessian))
}

