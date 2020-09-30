## ok, now let's try bsts, which has heavytailed options
## we use with the United States
library(stats)
library(MARSS)
library(forecast)
library(datasets)
library(demography)
library(bsts)

get_kt <- function(CODE,
                   years = 1950:2013,
                   ages = 0:100,
                   my.adjust = "e0")
{
    mx.prefix <- "~/Documents/hmd/hmd_statistics/death_rates/Mx_1x1/"
    expos.prefix <- "~/Documents/hmd/hmd_statistics/exposures/Exposures_1x1/"

    mort <- read.demogdata(paste0(mx.prefix, CODE, ".Mx_1x1.txt"),
               paste0(expos.prefix, CODE, ".Exposures_1x1.txt"),
               type = "mortality",
               label = CODE)
    min.mort.year <- min(years)
    max.mort.year <- max(years)
    if (min(mort$year) > 1950)
        min.mort.year <- min(mort$year)
    if (max(mort$year) < 2013)
        max.mort.year <- max(mort$year)
    years <- min.mort.year:max.mort.year

    mort.years <- extract.years(mort, years = years)
    mort.years.ages <- extract.ages(mort.years, ages = ages)
    lc.out <- lca(mort.years.ages, series = "total", adjust= my.adjust,
                  interpolate = TRUE)
    kt <- lc.out$kt
    return(kt)
}

usa_kt = as.vector(get_kt("USA"))
swe_kt = as.vector(get_kt("SWE"))


## semilocal
ss = list()
ss = AddSemilocalLinearTrend(ss, usa_kt)
semi_model = bsts(usa_kt,
             state.specification = ss,
             niter = 2000)
semi_pred = predict(semi_model, horizon = 50)
## local
ss = list()
ss = AddLocalLinearTrend(ss, usa_kt)
model = bsts(usa_kt,
             state.specification = ss,
             niter = 2000)
pred = predict(model, horizon = 50)
## student
ss = list()
ss = AddStudentLocalLinearTrend(ss, usa_kt)
stud_model = bsts(usa_kt,
             state.specification = ss,
             niter = 2000)
stud_pred = predict(stud_model, horizon = 50)



par(mfrow = c(3,2))
plot(model)
plot(pred, ylim = c(-500, 500))
plot(semi_model)
plot(semi_pred, ylim = c(-500, 500))
plot(stud_model)
plot(stud_pred, ylim = c(-500, 500))

sts = StructTS(usa_kt, type = "trend")

## let's compare the trends mu_t


burnin = 500
trend_model = colMeans(model$state.contributions[-(1:burnin),1,])
trend_semi_model = colMeans(semi_model$state.contributions[-(1:burnin),1,])
trend_stud_model = colMeans(stud_model$state.contributions[-(1:burnin),1,])

par(mfrow = c(1,1))
plot(trend_model, type = "l", col = "orange")
lines(trend_semi_model, col = "red")
lines(trend_stud_model, col = "blue")


plot(fitted(sts)[,"level"])
lines(trend_model, type = "l", col = "orange")

get_latent_state <- function(kt)
{
    ## local
    ss = list()
    ss = AddLocalLinearTrend(ss, kt)
    model = bsts(kt,
                 state.specification = ss,
                 niter = 2000)
    burnin = 500
    trend_model = colMeans(model$state.contributions[-(1:burnin),1,])
    return(trend_model)
    }


########################### let's see if we can do get bsts and StructTS doing same things ...

## let's do local linear trend model (this
y = swe_kt ## usa_kt
ss <- AddLocalLinearTrend(list(), y)
llt = bsts(y, niter = 2000, state.specification = ss)

## how do I get coefs?
## how do I get state?
burnin = 500
llt_state = colMeans(llt$state.contributions[-(1:burnin),1,])

## now let's compare to StructTS
sts_trend = StructTS(y, type = "trend")
sts_level = StructTS(y, type = "level")


plot(fitted(sts_trend)[,"level"], type = 'l')
lines(fitted(sts_level)[,"level"], type = 'l', lty = 2)
lines(llt_state, type = "l", col = "orange")
lines(y, col = "black", type = 'p')

## so bsts is way smoother than StructTS. Let's try DLM

############## let's try KFAS

## install.packages("KFAS")

library(KFAS)

y = usa_kt

model_gauss = SSModel(y ~ -1 +
                          SSMcustom
