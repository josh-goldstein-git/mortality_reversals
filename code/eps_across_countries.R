## we look at epsilon across countries

## we do marss lee-carter for our 9 countries and show difference in kt
## install.packages("demography")
library(demography)
library(data.table)

## countries (chosen manually)
files = system("ls ~/Documents/hmd/hmd_statistics/death_rates/Mx_1x1/*", intern = T)
files = gsub("^.*rates/Mx_1x1/", "", files)
files = gsub("\\..*$", "", files)
code.vec = files

code.vec = c("AUS", "AUT", "BEL", "CAN", "CZE", ## "DNK", ## "EST",
             "FIN", "FRATNP",
             "DEUTW", "ESP", "GBRTENW", "GRC", "HUN", "IRL", ## "ISL",
             "ISR", "ITA",
             "JPN", "NLD", "NOR", "NZL_NP", "POL", "RUS", "SWE", "TWN", "USA")
get_kt <- function(CODE,
                   years = 1950:2013, ## note this are bounds, works if less
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
    ## redefine years if not all are available
    if (min(mort$year) > 1950)
        min.mort.year <- min(mort$year)
    if (max(mort$year) < 2013)
        max.mort.year <- max(mort$year)
    years <- min.mort.year:max.mort.year
    mort.years <- extract.years(mort, years = years)
    mort.years.ages <- extract.ages(mort.years, ages = ages)
    ## do LC
    lc.out <- lca(mort.years.ages, series = "total", adjust= my.adjust,
                  interpolate = TRUE)
    kt <- lc.out$kt
    ## clean up and assign names
    kt.vec = as.vector(kt)
    names(kt.vec) <- time(kt)
    drift = (kt[length(kt)] - kt[1])/(length(kt)-1)
    eps = c(NA, diff(kt.vec) - drift)
    names(eps) <- time(kt)
    return(list(kt.vec = kt.vec,
           drift = drift,
           eps = eps))
}

year.vec = 1950:2018
eps.mat = matrix(NA, nrow = length(year.vec), ncol = length(code.vec))
colnames(eps.mat) <- code.vec
rownames(eps.mat) <- year.vec
for (i in 1:length(code.vec))
{
    my_code = code.vec[i]
    print(my_code)
    tmp = get_kt(my_code, years = year.vec)
    print(tmp$drift)
    eps.mat[names(tmp$eps),i] <- tmp$eps
}

## see which countries have eps


## now do analyis of eps.mat

R = cor(eps.mat[-1,], use = "pairwise")
round(R,1)

## let's use PCA to look at shocks

out = princomp(eps.mat[-1,], cor = TRUE, na.action = na.exclude)

