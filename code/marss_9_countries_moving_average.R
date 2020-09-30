
## we use moving average instead of MARSS

## (1) select countries and compute kt for them (should we use common bx, Q for Ron)
## (2) plot kt, smooths, and choose a moving average width (e.g., 5 years)
## (3) PCA of residuals with plots

## we do marss lee-carter for our 9 countries and show difference in kt
## install.packages("demography")
library(demography)
library(data.table)
library(MARSS)

## countries (chosen manually)
files = system("ls ~/Documents/hmd/hmd_statistics/death_rates/Mx_1x1/*", intern = T)
files = gsub("^.*rates/Mx_1x1/", "", files)
files = gsub("\\..*$", "", files)
code.vec = files

code.vec = c("AUS", "AUT", "BEL", "CAN", "CZE", "DNK", "EST", "FIN", "FRATNP",
             "DEUTW", "ESP", "GBRTENW", "GRC", "HUN", "IRL", ## "ISL",
             "ISR", "ITA",
             "JPN", "NLD", "NOR", "NZL_NP", "POL", "RUS", "SWE", "TWN", "USA")
## let's see the
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
    return(list(kt.vec = kt.vec,
           drift = drift,
           eps = eps))
}




## get kt on all countries
year.vec = 1950:2013
kt.mat = matrix(NA, nrow = length(year.vec), ncol = length(code.vec))
colnames(kt.mat) <- code.vec
rownames(kt.mat) <- year.vec
eps.mat = kt.mat
for (i in 1:length(code.vec))
{
    my_code = code.vec[i]
    print(my_code)
    tmp = get_kt(my_code, years = year.vec)
    kt.mat[names(tmp$kt.vec),i] <- tmp$kt.vec
    eps.mat[names(tmp$kt.vec),i] <- tmp$eps
}

write.table(x = kt.mat, "kt.mat")
foo = as.matrix(read.table("kt.mat"))
all.equal(kt.mat, foo)

kt.mat = as.matrix(read.table("kt.mat"))

## now get residuals of 5 year MA
n <- 7
kt_smu.mat <- apply(kt.mat, 2, filter, filter = rep(1/n, n))
kt_resid.mat <- kt.mat - kt_smu.mat
## let's plot the residuals
## opar <- par()
## get mean of residuals by year to get a sense of global shocks
global.vec = apply(kt_resid.mat, 1, mean, na.rm = T)

## eps.mat and kt_resid.mat, not super correlated
plot(eps.mat, kt_resid.mat)
cor(as.vector(eps.mat), as.vector(kt_resid.mat), use = "pairwise.complete.obs")
## [1] 0.5875669


## A = kt_resid.mat
A = eps.mat
par(mfrow = c(6,5), mar = c(0,0,1,0))
for (i in 1:ncol(A))
{
    plot(rownames(A), A[,i], main = colnames(A)[i], ylim = c(-10, 10), type  = "l")
    lines(rownames(A), global.vec, col = "grey", lwd = 1)
    lines(rownames(A), A[,i], main = colnames(A)[i])
    abline(v = 2003, col = "red") ## european heatwave
    abline(v = 1968+1, col = "blue") ## 1968 pandemix H3H2
}


cor(cbind(kt_resid.mat, global.vec), use = "pairwise.complete.obs")

## ok looks good

## get range of observations (not NA)
get.range = function(x, years)
{
    range(years[!is.na(x)])
}

t(apply(A, 2, get.range, years = rownames(A)))
##         [,1]   [,2]
## AUS     "1953" "2010"
## AUT     "1953" "2010"
## BEL     "1953" "2010"
## CAN     "1953" "2010"
## CZE     "1953" "2010"
## DNK     "1953" "2010"
## EST     "1962" "2010"
## FIN     "1953" "2010"
## FRATNP  "1953" "2010"
## DEUTW   "1959" "2010"
## ESP     "1953" "2010"
## GBRTENW "1953" "2010"
## GRC     "1984" "2010"
## HUN     "1953" "2010"
## IRL     "1953" "2010"
## ISL     "1953" "2010"
## ISR     "1986" "2010"
## ITA     "1953" "2010"
## JPN     "1953" "2010"
## NLD     "1953" "2010"
## NOR     "1953" "2010"
## NZL_NP  "1953" "2010"
## POL     "1961" "2010"
## RUS     "1962" "2010"
## SWE     "1953" "2010"
## TWN     "1973" "2010"
## USA     "1953" "2010"
######## note: the range is influenced by width of Moving Ave.
foo = t(apply(A, 2, get.range, years = rownames(A)))
## full.countries = rownames(foo[foo[,1] == "1953",])
start.year = as.numeric(foo[,1])
countries.1973 = rownames(foo[start.year <= 1973,])



####### pca

## let's only get the countries that don't have NAs
## B = A[,full.countries]
B = A[,countries.1973]
BB = B[!is.na(B[,1]),]

pca_result = prcomp(t(BB), scale = TRUE)

pc1 = pca_result$x[,1]
pc2 = pca_result$x[,2]
## let's plot PCA1 and PCA2
biplot(pca_result)

par(mfrow = c(1,1), mar = c(c(5, 4, 4, 2) + 0.1))
plot(pc1, pc2, type = 'n')
text(pc1, pc2, names(pc1))
## looks like nonsense to me
## part of it is pc2 is not so meaningufl
## but pc1 might be.
dotchart(sort(pc1))
## yes, this makes more sense.

## note, we're not getting the East West, North South that I was expecting.
## JPN very weird.
## English speaking countries very similar
## could drop HUN and ISL?



############## I like the correlation plot

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r <- abs(cor(x, y, use = "pairwise.complete.obs"))
         txt <- format(c(r, 0.123456789), digits = digits)[1]
         txt <- paste0(prefix, txt)
         if(missing(cex.cor))
             cex.cor <- 0.8/strwidth(txt)
         text(0.5, 0.5, txt, cex = cex.cor * sqrt(r) * 1.5)
     }





## rearrange

library(corrr)
R = correlate(BB, use = "pairwise")
## reR <- rearrange(R, absolute = TRUE)
reR <- rearrange(R, absolute = FALSE)
RR = as.matrix(reR)
mode(RR) <- "numeric"
round(RR,1)
my_order = colnames(RR)[-1] ## remove first one ("rowname")

pairs(BB[,my_order],##
##      [, my_order], ## lower.panel = panel.smooth,
      upper.panel = panel.cor,
      digits = 1,
           gap=0, row1attop=FALSE)

############# let's look at world before 1950-1970, 1970-1990, and 1990 to presetn





## try to reproduce order using PCA
getAnywhere("rearrange.cor_df")
out = seriation::seriate(BB, method = "PCA")
names(out[[2]])

## standardize to mean 0, sd = 1
myAs = apply(A, 2, scale)
## now get covariance
As.cov = cov(As, use='pairwise.complete.obs')
##
As.eigen = eigen(As.cov)

## look at eigen values
As.eigen$values[1:10]
 ## [1] 6.6556089 3.0009466 2.3221483 2.0775984 1.9039416 1.4864887 1.2973167
 ## [8] 1.1785570 0.9234722 0.8000095
## seems like 1st i sbig but after that gradual decline (maybe 2 is bigger than the rest)

phi = As.eigen$vec[,1:2]
phi = -phi
rownames(phi) = colnames(A)

## do some plots

## pdf("nt_pairs.pdf")

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r <- abs(cor(x, y, use = "pairwise.complete.obs"))
         txt <- format(c(r, 0.123456789), digits = digits)[1]
         txt <- paste0(prefix, txt)
         if(missing(cex.cor))
             cex.cor <- 0.8/strwidth(txt)
         text(0.5, 0.5, txt, cex = cex.cor * sqrt(r) * 1.5)
     }

## my_order = c("FRATNP", "ITA", "ESP", "SWE", "GBRTENW", "USA", "CAN", "JPN", "RUS", "DEUTW")
pairs(kt_resid.mat,##
##      [, my_order], ## lower.panel = panel.smooth,
      upper.panel = panel.cor,
           gap=0, row1attop=FALSE)
## dev.off()


foo = kt_resid.mat

## rearrange

library(corrr)
R = correlate(foo, use = "pairwise")
reR <- rearrange(R, absolute = TRUE)
RR = as.matrix(reR)
mode(RR) <- "numeric"
round(RR,1)
my_order = colnames(RR)[-1]

pairs(foo[,my_order],##
##      [, my_order], ## lower.panel = panel.smooth,
      upper.panel = panel.cor,
           gap=0, row1attop=FALSE)




###################################### end of moving ave analysis




fra_kt = as.vector(get_kt("FRATNP"))

fwrite(as.list(fra_kt), "../data/fra_kt.txt")

my_plot <- function(my_code)
{
    ## get kt
    kt = get_kt(my_code)
    year = as.vector(time(kt))
    kt = as.vector(kt)
    ## RWD
    mod.list.1 <- list(B = matrix(1),
                       U = matrix("d"),
                       Q = matrix("q"),
                       Z = matrix(1),
                       A = matrix(0),
                       R = matrix(0),
                       x0 = matrix("mu"),
                       tinitx = 0)
    lc.1 = MARSS(kt, model = mod.list.1)
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

    plot(year, as.vector(lc.1$states), type= "l", ylab = "kt", xlab = "year")
    lines(year, as.vector(lc.2$states), col = "red", lty =1, lwd = 1)
    legend("topright", legend = c("observed kt", "latent kt"),
           col = c("black", "red"), lty = c(1,1), lwd = 1:1)
    title(my_code)
}

my_plot_recent <- function(my_code)
{
    ## get kt
    kt = get_kt(my_code)
    year = as.vector(time(kt))
    kt = as.vector(kt)
    ## RWD
    mod.list.1 <- list(B = matrix(1),
                       U = matrix("d"),
                       Q = matrix("q"),
                       Z = matrix(1),
                       A = matrix(0),
                       R = matrix(0),
                       x0 = matrix("mu"),
                       tinitx = 0)
    lc.1 = MARSS(kt, model = mod.list.1)
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

    s = year >= 1990
    plot(year[s], as.vector(lc.1$states)[s], type= "l", ylab = "kt", xlab = "year")
    lines(year[s], as.vector(lc.2$states)[s], col = "red", lty =1, lwd = 1)
    legend("topright", legend = c("observed kt", "latent kt"),
           col = c("black", "red"), lty = c(1,1), lwd = 1:1)
    title(my_code)
}


my_plot_nt <- function(my_code)
{
    ## get kt
    kt = get_kt(my_code)
    year = as.vector(time(kt))
    kt = as.vector(kt)
    ## RWD
    mod.list.1 <- list(B = matrix(1),
                       U = matrix("d"),
                       Q = matrix("q"),
                       Z = matrix(1),
                       A = matrix(0),
                       R = matrix(0),
                       x0 = matrix("mu"),
                       tinitx = 0)
    lc.1 = MARSS(kt, model = mod.list.1)
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

    nt = as.vector(lc.2$states) - as.vector(lc.1$states)

    plot(nt, type= "l", ylim = c(-5, 5), ylim = )
    title(my_code)
}


pdf("kt_panel_fig.pdf", height = 10, width = 10)
par(mfrow = c(3,3))
for (i in 1:length(code.vec))
{
    my_code = code.vec[i]
    print(my_code)
    ##     my_plot_reveal(my_code)
        my_plot(my_code)
}
dev.off()
system('open kt_panel_fig.pdf')

## Let's zoom in on last 20 years


pdf("kt_resid_panel_fig.pdf", height = 10, width = 10)
par(mfrow = c(3,3))
for (i in 1:length(code.vec))
{
    my_code = code.vec[i]
        my_plot_reveal(my_code)
   ##     my_plot(my_code)
}
dev.off()

pdf("kt_recent_panel_fig.pdf", height = 10, width = 10)
par(mfrow = c(3,3))
for (i in 1:length(code.vec))
{
    my_code = code.vec[i]
    my_plot_recent(my_code)
   ##     my_plot(my_code)
}
dev.off()


## looks great -- although some did not converge

## let's do an AIC tabel


my_model <- function(my_code)
{
    ## get kt
    kt = get_kt(my_code)
    kt = as.vector(kt)
    ## RWD
    mod.list.1 <- list(B = matrix(1),
                       U = matrix("d"),
                       Q = matrix("q"),
                       Z = matrix(1),
                       A = matrix(0),
                       R = matrix(0),
                       x0 = matrix("mu"),
                       tinitx = 0)
    lc.1 = MARSS(kt, model = mod.list.1)
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
    return(list(lc.1 = lc.1, lc.2 = lc.2))
}

## my_list = vector("list", length(code.vec))
## names(my_list) = code.vec
## for (i in 1:length(code.vec))
## {
##     my_code = code.vec[i]
##     my_list[i] = my_model(my_code)
##    ##     my_plot(my_code)
## }

## lapply(my_list[1], AIC)


my_plot_nt("FRATNP")

my_plot_nt <- function(my_code)
{
    ## get kt
    kt = get_kt(my_code)
    year = as.vector(time(kt))
    kt = as.vector(kt)
    ## RWD
    mod.list.1 <- list(B = matrix(1),
                       U = matrix("d"),
                       Q = matrix("q"),
                       Z = matrix(1),
                       A = matrix(0),
                       R = matrix(0),
                       x0 = matrix("mu"),
                       tinitx = 0)
    lc.1 = MARSS(kt, model = mod.list.1)
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

    nt = as.vector(lc.2$states) - as.vector(lc.1$states)

    plot(year, nt, type= "l", ylim = c(-5, 5))
    title(my_code)
    return(nt)
}

par(mfrow = c(2,2))
nt_fra = my_plot_nt("FRATNP")
abline(v = 2004, col = 'grey')
nt_ita = nt_ita = my_plot_nt("ITA")
abline(v = 2004, col = 'grey')
nt_esp = my_plot_nt("ESP")
abline(v = 2004, col = 'grey')
plot(nt_fra, nt_ita)
points(nt_fra, nt_esp, col = "grey", pch = 19)

pairs(cbind(nt_fra, nt_ita, nt_esp))

abline(v = 2004, col = 'grey')


nt_can = my_plot_nt("CAN")
nt_usa = my_plot_nt("USA")
## nt_germany = my_plot_nt("DEUTW")
nt_swe = my_plot_nt("SWE")
nt_uk = my_plot_nt("GBRTENW")
nt_jpn = my_plot_nt("JPN")

nt.mat = cbind(nt_fra, nt_ita, nt_esp, nt_swe, nt_uk, nt_usa, nt_can, nt_jpn)
pairs(nt.mat)

cor.mat = cor(nt.mat)
round(cor.mat, 1)


panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r <- abs(cor(x, y))
         txt <- format(c(r, 0.123456789), digits = digits)[1]
         txt <- paste0(prefix, txt)
         if(missing(cex.cor))
             cex.cor <- 0.8/strwidth(txt)
         text(0.5, 0.5, txt, cex = cex.cor * r)
     }


pdf("nt_pairs.pdf")
pairs(nt.mat, ## lower.panel = panel.smooth,
      upper.panel = panel.cor,
           gap=0, row1attop=FALSE)
dev.off()
