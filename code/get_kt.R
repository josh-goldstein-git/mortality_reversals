library(demography)
library(data.table)
code.vec = c("USA", "JPN", "RUS", "ITA", "FRATNP", "SWE", "ESP", "GBRTENW", "DEUTW")



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

kt = get_kt("JPN")
my.list = vector("list", length(code.vec))
names(my.list) = code.vec
for (i in 1:length(code.vec))
{
    print(code.vec[i])
    kt = get_kt(code.vec[i])
    my.list[[i]] = data.table(year = as.vector(time(kt)), kt = as.vector(kt))
}

dt = rbindlist(my.list, idcol = TRUE)

xt = dt[, xtabs(kt ~ year + .id)]

foo = as.data.frame.matrix(xt)

out = as.data.table(foo)

fwrite(out, "kt_dt.csv")

bar = fread("kt_dt.csv")

