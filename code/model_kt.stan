data {
  int n;
  vector[n] y;
}

parameters {
  real d;
  real<lower=0> sigma_eps;
  real<lower=0> sigma_obs;
  vector[n] x;
}


model {
  d ~ normal(0,5);
  sigma_eps ~ normal(0,1);
  sigma_obs ~ normal(0,1);
  y ~ normal(x, sigma_obs);
  
  for (i in 2:n) {
    x[i] ~ normal(x[i-1] + d, sigma_eps);
  }
  // no forecast for now  
}

