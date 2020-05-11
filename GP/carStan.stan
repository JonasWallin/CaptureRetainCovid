data {
  int<lower=1> t;
  int N[t];
}

parameters {
  vector[t]       theta;
  real tau;
}
model {
    // hierarchical priors
  for(i in 2:t) {
    // non-centered parameterization
    theta[i] ~ normal(theta[i-1], 1/tau);
  }
      N ~ poisson_log(theta);
}


