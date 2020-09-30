## let's see if we can fit the RWD model with stan

library(rstan)
library(MARSS)


kt = scan("../../data/fra_kt.txt", sep = ",")

## let's do rwd with MARSS
mod.list.1 <- list(B = matrix(1),
                   U = matrix("d"),
                   Q = matrix("q"),
                   Z = matrix(1),
                   A = matrix(0),
                   R = matrix(0),
                   x0 = matrix("mu"),
                   tinitx = 0)
lc.1 = MARSS(kt, model = mod.list.1)
##       Estimate
## U.d      -1.93
## Q.q       5.09
## x0.mu    60.83
## Initial states (x0) defined at t=0



standata = list(y = kt, n = length(kt))

fit <- stan(file = "rwd.stan", data = standata,
            warmup = 100, iter = 2000, chains = 4)

stopifnot(is.converged(fit))
##         mean se_mean   sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
## sigma   2.34    0.01 0.21   1.98   2.19   2.33   2.47   2.81  1509    1
## d      -1.93    0.01 0.30  -2.50  -2.14  -1.93  -1.71  -1.34  1168    1
## lp__  -83.44    0.03 0.97 -85.98 -83.82 -83.14 -82.73 -82.45  1017    1

## 2.34^2
## [1] 5.4756

## I think the difference is because x0 is free param in MARSS



