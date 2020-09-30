data {
  int n_years;
  int n_countries;
  matrix[n_years, n_countries] Y;
}

parameters {
  real d_vec[n_countries];
  real<lower=0> sigma_eps_vec[n_countries];
  real<lower=0> sigma_obs_vec[n_countries];
  matrix[n_years,n_countries] X;  
}


model {
  for (i in 1:n_countries) {
    d_vec[i] ~ normal(0,5);
    sigma_eps_vec[i] ~ normal(0,1);
    sigma_obs_vec[i] ~ normal(0,1);
  y ~ normal(x, sigma_obs);
  
  for (i in 2:n) {
    x[i] ~ normal(x[i-1] + d, sigma_eps);
  }
  // no forecast for now  
}

