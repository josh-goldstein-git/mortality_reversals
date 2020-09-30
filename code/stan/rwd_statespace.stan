data {
  int<lower=1> n;
  vector[n] y;
}

parameters {
  vector[n] mu;
  real<lower=0> sigma_eps; // SD of innovation epsilon_t (or q in MARSS)
  real<lower=0> sigma_obs; // SD of obs
  real d; // drift
}

model {
  mu[1] ~ normal(y[1], sigma_obs);
  for(t in 2:n) // random_walk on mu_t with SD of level
    mu[t] ~ normal(mu[t-1] + d, sigma_eps);

  y ~ normal(mu, sigma_obs);
  
}
