data {
  int<lower=1> n;
  vector[n] y;
}

parameters {
  real<lower=0> sigma; // SD of innovation epsilon_t (or q in MARSS)
  real d; // drift
}

model {
  for(t in 2:n) // random_walk on mu_t with SD of level
    y[t] ~ normal(y[t-1] + d, sigma);
}
