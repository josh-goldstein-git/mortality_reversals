library(stats)
library(MARSS)
library(forecast)
library(datasets)
library(demography)
library(bsts)

code.vec = c("USA", "JPN", "CAN", "ITA", "FRATNP", "SWE", "ESP", "GBRTENW")

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
    return(list(model = model, trend_model = trend_model))
}

model_latent <- function(kt)
{
    ## local
    ss = list()
    ss = AddLocalLinearTrend(ss, kt)
    model = bsts(kt,
                 state.specification = ss,
                 niter = 2000)
    burnin = 500
    trend_model = colMeans(model$state.contributions[-(1:burnin),1,])
    return(list(model = model))
}



bsts_my_plot <- function(my_code)
{
    ## get kt
    kt = get_kt(my_code)
    year = as.vector(time(kt))
    kt = as.vector(kt)
    ## RWD
    ## RWD + obs error "R"
    kt_latent = get_latent_state(kt)
    plot(year, kt, type= "l", ylab = "kt", xlab = "year")
    lines(year, kt_latent, col = "red", lty =1, lwd = 1)
    legend("topright", legend = c("observed kt", "latent kt"),
           col = c("black", "red"), lty = c(1,1), lwd = 1:1)
    title(my_code)
}


## let's do all of the countries as a list

out_list = vector("list", length(code.vec))
names(out_list) = code.vec
for (i in 1:length(code.vec))
{
    my_code = code.vec[i]
        print(my_code)
    kt = get_kt(my_code)
    year = as.vector(time(kt))
    kt = as.vector(kt)
    out_list[[i]] = model_latent(kt)
}

## now get the latent_states

get_trend = function(country)
{
        burnin = 500
        trend = colMeans(country$model$state.contributions[-(1:burnin),1,])
        return(trend)
}

trend_list = lapply(out_list, get_trend)
kt_list = vector("list", length(code.vec))
for (i in 1:length(code.vec))
{
    my_code = code.vec[i]
    print(my_code)
    my_kt = as.vector(get_kt(my_code))
    kt_list[[i]] = my_kt
}
names(kt_list) = names(trend_list)

nt_list = vector("list", length(code.vec))
names(nt_list) = names(trend_list)
for (i in 1:length(code.vec))
{
    this_obs_kt = kt_list[[i]]
    this_latent_kt = trend_list[[i]]
    this_nt = this_latent_kt - this_obs_kt
    nt_list[[i]] = this_nt
}


par(mfrow = c(3,3))
for (i in 1:length(code.vec))
{
    my_code = code.vec[i]
    my_kt = kt_list[[my_code]]
    my_latent_kt = trend_list[[my_code]]
    my_nt = nt_list[[my_code]]
    plot(my_kt, type = 'l')
    lines(my_latent_kt, col = "orange")
    title(my_code)
    }

foo = cbind("USA" = nt_list$USA,
            "CAN" = nt_list$CAN,
            "JPN" = nt_list$JPN,
            "ITA" = nt_list$ITA,
            "FRATNP" = nt_list$FRATNP,
            "ESP" = nt_list$ESP,
            "SWE" = nt_list$SWE)

pairs(foo)
