## let's see if we can fit the RWD model with stan

library(rstan)
library(MARSS)


kt = scan("../../data/fra_kt.txt", sep = ",")

## let's do rwd with MARSS
## RWD + obs error "R"
mod.list.2 <- list(B = matrix(1),
                   U = matrix("d"),
                   Q = matrix("q"),
                   Z = matrix(1),
                   A = matrix(0),
                   R = matrix("r"),
                   x0 = matrix("mu"),
                   tinitx = 0)
lc.2 = MARSS(kt, model = mod.list.2)
##       Estimate
## R.r       1.82
## U.d      -1.94
## Q.q       1.26
## x0.mu    61.45
## Initial states (x0) defined at t=0



standata = list(y = kt, n = length(kt))

fit <- stan(file = "rwd_statespace.stan", data = standata,
            warmup = 2000, iter = 4000, chains = 4)
## mu[63]    -60.84    0.01 0.88  -62.59 -61.43 -60.84 -60.23 -59.15  5777 1.00
## mu[64]    -62.76    0.01 1.02  -64.78 -63.43 -62.78 -62.09 -60.73  5709 1.00
## sigma_eps   1.30    0.01 0.28    0.85   1.10   1.27   1.46   1.93  1025 1.00
## sigma_obs   1.30    0.01 0.23    0.83   1.16   1.30   1.45   1.76  1201 1.00
## d          -1.94    0.00 0.17   -2.27  -2.05  -1.94  -1.83  -1.59  7292 1.00
## lp__      -93.87    0.26 8.96 -112.39 -99.52 -93.62 -87.71 -77.14  1180 1.01


mu = get_posterior_mean(fit, par = "mu")[, "mean-all chains"]
mu = as.vector(mu)

mu_marss = as.vector(lc.2$states)
plot(mu, type = 'l')
lines(mu_marss, col = 'red')
points(kt, cex = .5)
## bingo. Works fine



fit <- stan(file = "rwd_statespace_student_t.stan", data = standata,
            warmup = 4000, iter = 8000, chains = 4)

mu_student = get_posterior_mean(fit, par = "mu")[, "mean-all chains"]


plot(mu, type = 'l')
lines(mu_marss, col = 'red')
points(kt, cex = .5)
lines(mu_student, col = 'blue')
