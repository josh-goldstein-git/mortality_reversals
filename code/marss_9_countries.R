## we do marss lee-carter for our 9 countries and show difference in kt

library(data.table)
library(demography)

code.vec = c("USA", "JPN", "RUS", "ITA", "FRATNP", "SWE", "ESP", "GBRTENW", "DEUTW")

code.vec = c("AUS", "AUT", "BEL", "CAN", "CZE", ## "DNK", "EST", "FIN",
             "FRATNP",
             "DEUTW", "ESP", "GBRTENW", ## "GRC", "HUN", "IRL", ## "ISL",
             ## "ISR",
             "ITA",
             "JPN",
             "NLD", "NOR", "NZL_NP", "POL", "RUS", "SWE", "TWN", "USA")


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

fra_kt = as.vector(get_kt("FRATNP"))
library(data.table)
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
