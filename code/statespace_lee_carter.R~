## Let's do statespace LC with different forms of the model
## (1) deterministic drift with random walk innovations
## (2) deterministic drift with random walk innovations and observation error
## (3) random drift with random walk innovations and observation error

library(stats)
library(MARSS)
library(forecast)
library(datasets)
library(demography)

## Let's start with Sweden

mx.prefix <- "~/Documents/hmd/hmd_statistics/death_rates/Mx_1x1/"
expos.prefix <- "~/Documents/hmd/hmd_statistics/exposures/Exposures_1x1/"
swe.mort <- read.demogdata(paste0(mx.prefix, "SWE.Mx_1x1.txt"),
               paste0(expos.prefix, "SWE.Exposures_1x1.txt"),
               type = "mortality",
               label = "SWE")
mort <- swe.mort

years = 1950:2013
min.mort.year <- min(years)
max.mort.year <- max(years)
if (min(mort$year) > 1950)
    min.mort.year <- min(mort$year)
if (max(mort$year) < 2013)
    max.mort.year <- max(mort$year)
years <- min.mort.year:max.mort.year
ages = 0:100
my.adjust = "e0"
mort.years <- extract.years(mort, years = years)
mort.years.ages <- extract.ages(mort.years, ages = ages)
lc.out <- lca(mort.years.ages, series = "total", adjust= my.adjust,
              interpolate = TRUE)
    kt <- lc.out$kt
    n <- length(kt)
    drift <- -(kt[1] - kt[n]) / (n-1)
    eps <- diff(kt) - drift
    eps.std <- eps/abs(drift)
    eps.std <- unclass(eps.std)

kt <- as.vector(kt)
swe_kt = kt

## Now let's do statespace models
## cannonical form : x[t] = B * x[t-1] + u + wt ; wt ~ N(0,Q)
## cannonical form : y[t] = Z * x[t]  + a + vt; vt ~ N(0, R)

## model 1: LC RWD

mod.list.1 <- list(B = matrix(1),
                 U = matrix("d"),
                 Q = matrix("q"),
                 Z = matrix(1),
                 A = matrix(0),
                 R = matrix(0),
                 x0 = matrix("mu"),
                 tinitx = 0)
lc.1 = MARSS(kt, model = mod.list.1)


## model 2: LC RWD + observation errro
mod.list.2 <- list(B = matrix(1),
                 U = matrix("d"),
                 Q = matrix("q"),
                 Z = matrix(1),
                 A = matrix(0),
                 R = matrix("r"),
                 x0 = matrix("mu"),
                 tinitx = 0)
lc.2 = MARSS(kt, model = mod.list.2)

## model 3: LC RWD + observation error + drift error (use structTS)
lc.3 = StructTS(x = kt, type = "trend")
## Variances:
##   level    slope  epsilon
## 0.51319  0.01591  1.48278

AIC(lc.1) ## [1] 265.9943
AIC(lc.2) ## [1] 258.5834 ## better

MARSSparamCIs(lc.1)
##       ML.Est Std.Err low.CI up.CI
## U.d    -1.73   0.232  -2.19 -1.28
## Q.q     3.40   0.601   2.22  4.58
## x0.mu  51.56   1.859  47.91 55.20
## Initial states (x0) defined at t=0

MARSSparamCIs(lc.2)
##       ML.Est Std.Err low.CI up.CI
## R.r     1.04   0.392  0.267  1.81
## U.d    -1.73   0.142 -2.003 -1.45
## Q.q     1.24   0.498  0.267  2.22
## x0.mu  51.33   1.401 48.579 54.07
## Initial states (x0) defined at t=0

plot(as.vector(lc.1$states), type= "l")
lines(as.vector(lc.2$states), col = "red")
## way smoother

## but in this case ends on same thing.

## let's do USA thru 2017


mx.prefix <- "~/Documents/hmd/hmd_statistics/death_rates/Mx_1x1/"
expos.prefix <- "~/Documents/hmd/hmd_statistics/exposures/Exposures_1x1/"
usa.mort <- read.demogdata(paste0(mx.prefix, "USA.Mx_1x1.txt"),
               paste0(expos.prefix, "USA.Exposures_1x1.txt"),
               type = "mortality",
               label = "USA")
mort <- usa.mort
## years = 1950:2017
years = 1933:2017
min.mort.year <- min(years)
max.mort.year <- max(years)
if (min(mort$year) > 1950)
    min.mort.year <- min(mort$year)
if (max(mort$year) < 2013)
    max.mort.year <- max(mort$year)
years <- min.mort.year:max.mort.year
ages = 0:100
my.adjust = "e0"
mort.years <- extract.years(mort, years = years)
mort.years.ages <- extract.ages(mort.years, ages = ages)
lc.out <- lca(mort.years.ages, series = "total", adjust= my.adjust,
              interpolate = TRUE)
    kt <- lc.out$kt
    n <- length(kt)
    drift <- -(kt[1] - kt[n]) / (n-1)
    eps <- diff(kt) - drift
    eps.std <- eps/abs(drift)
    eps.std <- unclass(eps.std)
kt <- as.vector(kt)
usa_kt = kt

## model 1: LC RWD
lc.1 = MARSS(kt, model = mod.list.1)
## model 2: LC RWD + observation errro
## lc.2 = MARSS(kt, model = mod.list.2, control=list(maxit=10000))

lc.2 = MARSS(kt[1:60], model = mod.list.2)

plot(as.vector(lc.1$states), type= "l")
lines(as.vector(lc.2$states), col = "red")

lc.3 = StructTS(x = kt, type = "trend")
tsdiag(lc.3)

filt = KalmanRun(kt, lc.3$model)
plot(kt, type = "l")
lines(filt$states[,1], lty = 2, col = "red")

## plot state of drift
par(mfrow = c(2,1))
plot(filt$states[,1], lty = 2, col = "red")
plot(filt$states[,2], lty = 2, col = "red")
## I'm confused, why is state[,2] so variable?


########### let's revist SWE
n = length(swe_kt)
swe_lc.1 = MARSS(swe_kt[-n], model = mod.list.1)
swe_lc.2 = MARSS(swe_kt[-n], model = mod.list.2)
swe_lc.3 = StructTS(x = swe_kt[-n], type = "trend")
tsdiag(lc.3)
par(mfrow = c(1,1))
plot(as.vector(swe_lc.1$states), type= "l")
lines(as.vector(swe_lc.2$states), col = "red")

## we clearly get a different take-off here.

## now revisit USA and see if we can get a model 2 to work with different window

my_usa_kt = usa_kt[years %in% 1946:2014]
usa_lc.1 = MARSS(my_usa_kt, model = mod.list.1)
usa_lc.2 = MARSS(my_usa_kt, model = mod.list.2)
##
par(mfrow = c(1,1))
plot(as.vector(usa_lc.1$states), type= "l")
lines(as.vector(usa_lc.2$states), col = "red")

## I think that there's not enough negative auto-correlation in the US
## data to make a difference here.
## 1. -- Can check with some other countries
## 2. -- (advanced) can vary the form of the shock in Stan to take into account
## the nature of the shocks (fat tailed, one-sided etc)

## let's do FRANCE

mx.prefix <- "~/Documents/hmd/hmd_statistics/death_rates/Mx_1x1/"
expos.prefix <- "~/Documents/hmd/hmd_statistics/exposures/Exposures_1x1/"
fratnp.mort <- read.demogdata(paste0(mx.prefix, "FRATNP.Mx_1x1.txt"),
               paste0(expos.prefix, "FRATNP.Exposures_1x1.txt"),
               type = "mortality",
               label = "FRATNP")
mort <- fratnp.mort
years = 1950:2017
## years = 1933:2017
min.mort.year <- min(years)
max.mort.year <- max(years)
if (min(mort$year) > 1950)
    min.mort.year <- min(mort$year)
if (max(mort$year) < 2013)
    max.mort.year <- max(mort$year)
years <- min.mort.year:max.mort.year
ages = 0:100
my.adjust = "e0"
mort.years <- extract.years(mort, years = years)
mort.years.ages <- extract.ages(mort.years, ages = ages)
lc.out <- lca(mort.years.ages, series = "total", adjust= my.adjust,
              interpolate = TRUE)
    kt <- lc.out$kt
    n <- length(kt)
    drift <- -(kt[1] - kt[n]) / (n-1)
    eps <- diff(kt) - drift
    eps.std <- eps/abs(drift)
    eps.std <- unclass(eps.std)
kt <- as.vector(kt)
fra_kt = kt

my_fra_kt = fra_kt ## [years %in% 1946:2014]
fra_lc.1 = MARSS(my_fra_kt, model = mod.list.1)
fra_lc.2 = MARSS(my_fra_kt, model = mod.list.2)
##
par(mfrow = c(1,1))
plot(as.vector(fra_lc.1$states), type= "l")
lines(as.vector(fra_lc.2$states), col = "red")


## create a function to do this and show for all 9 countries ... along with AIC
## this could highlight exceptional nature of US


## also to do is forecast -- hyndman funcitons show how to forecast StructTS
## could use to figure out how to do MARSS forecast (with level model)

