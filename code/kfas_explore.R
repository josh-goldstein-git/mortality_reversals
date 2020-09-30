## let's do random walk model

data("alcohol", package = "KFAS")
deaths = window(alcohol[, 2], end = 2007)
population = window(alcohol[, 6], end = 2007)
Zt = rbind(c(1,0))
Ht = matrix(NA)
Tt = matrix(c(1,0,1,1), 2,2)
Rt = cbind(c(1,0))
Qt = matrix(NA)
a1 = matrix(c(1,0), 2, 1)
P1 = matrix(0, 2,2)
P1inf = diag(2)

rate = deaths / population
model_gaussian = SSModel( rate ~ -1 +
                              SSMcustom(Z = Zt,
                                        T = Tt,
                                        R = Rt,
                                        Q = Qt,
                                        a1 = a1,
                                        P1 = P1,
                                        P1inf = P1inf),
                         H = Ht)

fit_gaussian = fitSSM(model_gaussian, inits = c(0, 0), method = "BFGS")
out_gaussian = KFS(fit_gaussian$model)
