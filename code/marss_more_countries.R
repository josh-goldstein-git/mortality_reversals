## we do marss lee-carter for our 9 countries and show difference in kt

library(data.table)
library(demography)
library(MARSS)


code.vec = c("USA", "JPN", "RUS", "ITA", "FRATNP", "SWE", "ESP", "GBRTENW", "DEUTW")

code.vec = c("AUS", "AUT", "BEL", "CAN", "CZE", ## "DNK", "EST", "FIN",
             "FRATNP",
             "DEUTW", "ESP", "GBRTENW", ## "GRC", "HUN", "IRL", ## "ISL",
             ## "ISR",
             "ITA",
             "JPN",
             "NLD", "NOR",
             "NZL_NP",
             "POL", "RUS", "SWE", "TWN", "USA")


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
    kt.vec = as.vector(kt)
    names(kt.vec) <- time(kt)
    return(kt.vec)
}


my_plot <- function(my_code)
{
    ## get kt
    kt = get_kt(my_code)
    year = names(kt)
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

my_plot("FRATNP")
my_plot("NOR")
my_plot("JPN")
my_plot("USA")


my_plot_recent <- function(my_code)
{
    ## get kt
    kt = get_kt(my_code)
    year = names(kt)
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
    year = names(kt)
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
par(mfrow = c(4,5))
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
par(mfrow = c(4,5))
for (i in 1:length(code.vec))
{
    my_code = code.vec[i]
        my_plot_nt(my_code)
   ##     my_plot(my_code)
}
dev.off()
system('open kt_resid_panel_fig.pdf')

pdf("kt_recent_panel_fig.pdf", height = 10, width = 10)
par(mfrow = c(3,3))
for (i in 1:length(code.vec))
{
    my_code = code.vec[i]
    my_plot_recent(my_code)
   ##     my_plot(my_code)
}
dev.off()

################## do as list

my_model <- function(my_code)
{
    ## get kt using Lee-Carter
    kt = get_kt(my_code)
    year = names(kt)
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
    return(list(year = year, kt = kt, lc.1 = lc.1, lc.2 = lc.2))
}


code.vec = c("FRATNP", "USA")
my_list = vector("list", length(code.vec))
names(my_list) = code.vec

for (i in 1:length(code.vec))
    my_list[[i]] = my_model(code.vec[i])


## now plot
my_plot_kt <- function(my_code, zoom_years = NULL, newtitle = NULL, ...)
{
    my_result = my_list[[my_code]]
    year = my_result$year
    lc.1 = my_result$lc.1
    lc.2 = my_result$lc.2
    ## now subset to zoom years, doing all if no zoom
    if (is.null(zoom_years))
        s = rep(TRUE, length(year))
    if (!is.null(zoom_years))
        s = year %in% zoom_years
    plot(year[s], as.vector(lc.1$states)[s], type= "l", ylab = "kt", xlab = "year", ...)
    lines(year[s], as.vector(lc.2$states)[s], col = "red", lty =1, lwd = 1)
    legend("topright", legend = c("observed kt", "latent kt"),
           col = c("black", "red"), lty = c(1,1), lwd = 1:1, bty = "n")
    if(is.null(newtitle))
        title(my_code)
    if(!is.null(newtitle))
        title(newtitle)
}

my_plot_nt <- function(my_code, zoom_years = NULL, newtitle = NULL, ...)
{
    my_result = my_list[[my_code]]
    year = my_result$year
    lc.1 = my_result$lc.1
    lc.2 = my_result$lc.2
    ## now subset to zoom years, doing all if no zoom
    if (is.null(zoom_years))
        s = rep(TRUE, length(year))
    if (!is.null(zoom_years))
        s = year %in% zoom_years


    nt = as.vector(lc.1$states) - as.vector(lc.2$states)
    names(nt) = year

    plot(year[s], nt[s], type= "l", ...)
    if(is.null(newtitle))
        title(my_code)
    if(!is.null(newtitle))
        title(newtitle)
    return(nt)
}


code.vec = c("AUS", "AUT", "BEL", "CAN", "CZE", ## "DNK", "EST", "FIN",
             "FRATNP",
             "DEUTW", "ESP", "GBRTENW", ## "GRC", "HUN", "IRL", ## "ISL",
             ## "ISR",
             "ITA",
             "JPN",
             "NLD", "NOR",
             "NZL_NP",
             "POL", "RUS", "SWE", "TWN", "USA")
my_list = vector("list", length(code.vec))
names(my_list) = code.vec

for (i in 1:length(code.vec))
    my_list[[i]] = my_model(code.vec[i])


for (i in 1:length(code.vec))
{
    my_plot_nt(code.vec[i], zoom_years = 1965:1985, ylim = c(-5, 5))
}


for (i in 1:length(code.vec))
{
    my_plot_kt(code.vec[i], zoom_years = 1990:2005)
}


## figure 1, France
pdf("france_example.pdf", width = 6, height = 6)
par(mfrow = c(1,1))
my_plot_kt("FRATNP", newtitle = "Lee-Carter mortality time trend (kt) in France")
abline(v = c(1968, 2003), lty = 2)
text(x = c(1968, 2003), y = c(-50, -50), c("1968\n Flu", "2003\n Heatwave"), pos = 2)
dev.off()
system("open france_example.pdf")

## ok, let's do a panel plot
pdf("kt_panel_plot.pdf", width = 10, height = 12)
par(mfrow = c(4,5))
for (i in 1:length(code.vec))
{
    my_plot_kt(code.vec[i], xlim = c(1950, 2015))
}
dev.off()
system("open kt_panel_plot.pdf")

pdf("nt_panel_plot.pdf", width = 10, height = 12)
nt.mat = matrix(NA, nrow = length(1950:2015), ncol = length(code.vec))
rownames(nt.mat) = 1950:2015
colnames(nt.mat) = code.vec
par(mfrow = c(4,5))
for (i in 1:length(code.vec))
{
    nt = my_plot_nt(code.vec[i], zoom_years = NULL, ylim = c(-5, 5),
                    xlim = c(1950, 2015))
    if(code.vec[i] == "FRATNP")
    {
        abline(v = c(1968, 2003), lty = 2)
        text(x = c(1968, 2003), y = c(-50, -50), c("1968\n Flu", "2003\n Heatwave"), pos = 2)
        }
    nt.mat[names(nt),i] <- nt
}
dev.off()
system("open nt_panel_plot.pdf")


## now pairs plot

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r <- abs(cor(x, y, use = "pairwise"))
         txt <- format(c(r, 0.123456789), digits = digits)[1]
         txt <- paste0(prefix, txt)
         if(missing(cex.cor))
             cex.cor <- 0.8/strwidth(txt)
         text(0.5, 0.5, txt, cex = cex.cor * sqrt(r)*1.3)
     }


R = cor(nt.mat, use = "pairwise")
## diag(R) = NA
e = eigen(R)
o = order(-e$vec[,1])


## pdf("nt_pairs.pdf")
## colnames(nt.mat)
##  [1] "AUS"     "AUT"     "BEL"     "CAN"     "CZE"     "FRATNP"  "DEUTW"
##  [8] "ESP"     "GBRTENW" "ITA"     "JPN"     "NLD"     "NOR"     "NZL_NP"
## [15] "POL"     "RUS"     "SWE"     "TWN"     "USA"
nice_names = c("Austra\n-lia", "Austria", "Belgium", "Canada", "Czech", "France", "West\n Germany",
               "Spain", "U.K.", "Italy", "Japan", "Netherl.", "Norway", "New\n Zealand",
               "Poland", "Russia", "Sweden", "Taiwan", "U.S.")

pdf("nt_corr_plot.pdf", width = 12, height = 12)
pairs(nt.mat[, o],
      upper.panel = panel.cor,
      digits = 1,
      labels = nice_names[o],
      cex.labels = 1.2,
      gap=0, row1attop=FALSE)
dev.off()
system("open nt_corr_plot.pdf")



## looks great

## let's do an AIC tabel




my_model("FRATNP")

my_list = vector("list", length(code.vec))
names(my_list) = code.vec
for (i in 1:length(code.vec))
{
    my_code = code.vec[i]
    print(my_code)
    my_list[[i]] = my_model(my_code)
    ##    ##     my_plot(my_code)
}


lapply(my_list, AIC)

## now plot going through the list

plot.kt <- function(my.code)
{
    X <- my_list[[my.code]]
    year
    plot(year, as.vector(X$lc.1$states), type= "l", ylab = "kt", xlab = "year")
    lines(year, as.vector(X$lc.2$states), col = "red", lty =1, lwd = 1)
    legend("topright", legend = c("observed kt", "latent kt"),
           col = c("black", "red"), lty = c(1,1), lwd = 1:1)
    title(my.code)
}


plot.kt("FRATNP")

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
