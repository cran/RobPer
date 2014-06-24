### R code from vignette source 'RobPer_Vignette.Rnw'

###################################################
### code chunk number 1: RobPer_Vignette.Rnw:143-149
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
library(RobPer)
## NOTE:
## Calculations and figures related to real data examples are not
## available as code chunks due to copyright uncertainties.
## This concerns Sections 5.2 and 5.3 as well as Figures 1 and 7 to 10.


###################################################
### code chunk number 2: RobPer_Vignette.Rnw:469-486
###################################################
# Calculations for Figure 2
set.seed(923)
zr <- tsgen(ttype="sine", npoints=500, ytype="const", pf=150,
    ncycles=71, ps=4, redpart=0.2, alpha=0,SNR=2, interval=FALSE)
PP <- RobPer(zr, model="sine", regression="L2", weighting=TRUE,
    periods=1:100)
shapes <- betaCvMfit(PP)
myf <- function(x) dbeta(x, shapes[1], shapes[2])
stdf <- function(x) dbeta(x, 1,497/2)

# Panel 2(a)
par(mar=c(4,4,0.1,0.1))
plot(PP, xlab="Trial period", ylab="Periodogram", type="l", axes=FALSE)
axis(1)
axis(2, at=seq(0,0.06, by=0.02))
abline(h=qbeta(0.95^0.01, 1,497/2), lty=2)
box()


###################################################
### code chunk number 3: RobPer_Vignette.Rnw:490-499
###################################################
# Panel 2(b)
par(mar=c(4,4,0.1,0.1))
hist(PP, freq=FALSE, main=" ", xlab="Periodogram",ylab="Density",
    col="grey80", breaks=15, axes=FALSE)
axis(1, at=seq(0,0.06, by=0.02))
axis(2)
box()
curve(myf, add=TRUE)
curve(stdf, add=TRUE, lty=2)


###################################################
### code chunk number 4: RobPer_Vignette.Rnw:547-573
###################################################
# Figure 3
par(mar=c(3,3,0.3,0.1), mgp=c(2,1,0))
set.seed(12)
PP <- c(rbeta(45, shape1=4, shape2=15), runif(5, min=0.8, max=1))
hist(PP, freq=FALSE, breaks=30, ylim=c(0,7), col="grey90",
    main="", xlab="Periodogram bar")
# true parameters:
myf.true <- function(x) dbeta(x, shape1=4, shape2=15)
curve(myf.true, add=TRUE, lwd=2)
# method of moments:
par.mom <- betaCvMfit(PP, rob=FALSE, CvM=FALSE)
myf.mom <- function(x) dbeta(x, shape1=par.mom[1], shape2=par.mom[2])
curve(myf.mom, add=TRUE, lwd=2, lty=2)
# robust method of moments
par.rob <- betaCvMfit(PP, rob=TRUE, CvM=FALSE)
myf.rob <- function(x) dbeta(x, shape1=par.rob[1], shape2=par.rob[2])
curve(myf.rob, add=TRUE, lwd=2, lty=3)
# CvM distance minimization
par.CvM <- betaCvMfit(PP, rob=TRUE, CvM=TRUE)
myf.CvM <- function(x) dbeta(x, shape1=par.CvM[1], shape2=par.CvM[2])
curve(myf.CvM, add=TRUE, lwd=2, col="grey50")
# Searching for outliers...
abline(v=qbeta((0.95)^(1/50), shape1=par.CvM[1], shape2=par.CvM[2]), col="grey30")
legend("topright", col=c("black", "grey50","black", "black"),lwd=2, lty=c(1,1,3,2),
    legend=c("True", "CvM", "Robust moments", "Moments"))
box()


###################################################
### code chunk number 5: RobPer_Vignette.Rnw:709-714
###################################################
library("RobPer")
set.seed(22)
lightcurve <- tsgen(ttype = "sine", ytype = "peak", pf = 7,
    redpart = 0.1, s.outlier.fraction = 0.1, interval = TRUE,
    npoints = 200, ncycles = 25, ps = 20, SNR = 3, alpha = 0)


###################################################
### code chunk number 6: RobPer_Vignette.Rnw:721-723
###################################################
set.seed(22)
tt <- sampler(ttype = "sine", npoints = 200, ncycles = 25, ps = 20)


###################################################
### code chunk number 7: RobPer_Vignette.Rnw:726-727
###################################################
yf <- signalgen(tt, ytype = "peak", pf = 7)


###################################################
### code chunk number 8: RobPer_Vignette.Rnw:730-733
###################################################
temp <- lc_noise(tt, sig = yf, SNR = 3, redpart = 0.1, alpha = 0)
y <- temp$y
s <- temp$s


###################################################
### code chunk number 9: RobPer_Vignette.Rnw:736-738
###################################################
temp <- disturber(tt, y, s, ps = 20, s.outlier.fraction = 0.1,
    interval = TRUE)


###################################################
### code chunk number 10: RobPer_Vignette.Rnw:741-742
###################################################
all(cbind(tt, temp$y, temp$s) == lightcurve)


###################################################
### code chunk number 11: RobPer_Vignette.Rnw:750-760
###################################################
# Panel 4(a)
par(mar=c(3,3,0.3,0.3), mgp=c(2,1,0))
plot(lightcurve[,1], lightcurve[,2], pch=16, xlab="t", ylab="y",
    main=" ", ylim=range(c(lightcurve[,2]+lightcurve[,3],lightcurve[,2]-lightcurve[,3])),
    col=rgb(0,0,0,alpha=0.5), cex=0.7, axes=FALSE)
axis(1, at=(0:5)*100, labels=c(0, "", 200, "", 400, ""))
axis(2)
rect(lightcurve[,1], lightcurve[,2]+lightcurve[,3], lightcurve[,1],
    lightcurve[,2]-lightcurve[,3], border=rgb(0,0,0,alpha=0.5))
box()


###################################################
### code chunk number 12: RobPer_Vignette.Rnw:764-771
###################################################
# Panel 4(b)
par(mar=c(3,3,0.3,0.3), mgp=c(2,1,0))
plot(lightcurve[,1]%%7, lightcurve[,2], pch=16, col=rgb(0,0,0,alpha=0.5),
    xlab="t modulo 7", ylab="y", main=" ",
    ylim=range(c(lightcurve[,2]+lightcurve[,3],lightcurve[,2]-lightcurve[,3])), cex=0.7)
rect(lightcurve[,1]%%7, lightcurve[,2]+lightcurve[,3], lightcurve[,1]%%7,
    lightcurve[,2]-lightcurve[,3], border=rgb(0,0,0,alpha=0.5))


###################################################
### code chunk number 13: RobPer_Vignette.Rnw:775-781
###################################################
# Panel 4 (c)
par(mar=c(3,3,0.3,0.3), mgp=c(2,1,0))
hist(lightcurve[,1]%%20, xlab="t modulo 20", col="grey", main=" ", freq=FALSE)
dsin <- function(tt) (sin(2*pi*tt/20)+1)/20
curve(dsin, add=TRUE, lwd=2)
box()


###################################################
### code chunk number 14: RobPer_Vignette.Rnw:790-792
###################################################
PP <- RobPer(lightcurve, model = "splines", regression = "huber",
    weighting = FALSE, var1 = FALSE, periods = 1:50)


###################################################
### code chunk number 15: RobPer_Vignette.Rnw:796-799
###################################################
betavalues <- betaCvMfit(PP)
crit.val <- qbeta((0.95)^(1 / 50), shape1 = betavalues[1],
    shape2 = betavalues[2])


###################################################
### code chunk number 16: RobPer_Vignette.Rnw:802-808
###################################################
hist(PP, breaks = 20, freq = FALSE, xlim = c(0, 0.08), col = "grey",
    main = "", xlab="Periodogram")
betafun <- function(x) dbeta(x, shape1 = betavalues[1],
    shape2 = betavalues[2])
curve(betafun, add = TRUE, lwd = 2)
abline(v = crit.val, lwd = 2)


###################################################
### code chunk number 17: RobPer_Vignette.Rnw:811-817
###################################################
par.mom <- betaCvMfit(PP, rob = FALSE, CvM = FALSE)
myf.mom <- function(x) dbeta(x, shape1 = par.mom[1], shape2 = par.mom[2])
curve(myf.mom, add = TRUE, lwd = 2, lty = 2)
crit.mom <- qbeta((0.95)^(1 / 50), shape1 = par.mom[1],
    shape2 = par.mom[2])
abline(v = crit.mom, lwd = 2, lty = 2)


###################################################
### code chunk number 18: RobPer_Vignette.Rnw:820-829
###################################################
par.rob <- betaCvMfit(PP, rob = TRUE, CvM = FALSE)
myf.rob <- function(x) dbeta(x, shape1 = par.rob[1], shape2 = par.rob[2])
curve(myf.rob, add = TRUE, lwd = 2, lty = 3)
crit.rob <- qbeta((0.95)^(1 / 50), shape1 = par.rob[1],
    shape2 = par.rob[2])
abline(v = crit.rob, lwd = 2, lty = 3)
legend("topright", lty = 1:3, legend = c("CvM", "Moments",
    "Robust moments"), bg = "white", lwd = 2)
box()


###################################################
### code chunk number 19: RobPer_Vignette.Rnw:832-837
###################################################
plot(1:50, PP, xlab = "Trial period", ylab = "Periodogram", main = "",
    type = "l")
abline(h = crit.val, lwd = 2)
text(7, PP[7]-0.002,7, pos=4)
text(14, PP[14]+0.002,14, pos=4)


###################################################
### code chunk number 20: RobPer_Vignette.Rnw:844-857
###################################################
# Panel 5(a)
par(mar=c(3,3,1,0.7), mgp=c(2,1,0))
hist(PP, breaks = 20, freq = FALSE, xlim = c(0, 0.08), col = "grey",
    main = "", xlab="Periodogram")
curve(betafun, add = TRUE, lwd = 2)
abline(v = crit.val, lwd = 2)
curve(myf.mom, add = TRUE, lwd = 2, lty = 2)
abline(v = crit.mom, lwd = 2, lty = 2)
curve(myf.rob, add = TRUE, lwd = 2, lty = 3)
abline(v = crit.rob, lwd = 2, lty = 3)
legend("topright", lty = 1:3, legend = c("CvM", "Moments", "Robust moments"),
    bg = "white", lwd = 2)
box()


###################################################
### code chunk number 21: RobPer_Vignette.Rnw:861-868
###################################################
# Panel 5(b)
par(mar=c(3,3,1,0.1), mgp=c(2,1,0))
plot(1:50, PP, xlab = "Trial period", ylab = "Periodogram", main = "",
    type = "l")
abline(h = crit.val, lwd = 2)
text(7, PP[7]-0.002,7, pos=4)
text(14, PP[14]+0.002,14, pos=4)


###################################################
### code chunk number 22: RobPer_Vignette.Rnw:877-879
###################################################
PP <- RobPer(lightcurve, model = "splines", regression = "L2",
    weighting = FALSE, var1 = FALSE, periods = 1:50)


###################################################
### code chunk number 23: RobPer_Vignette.Rnw:886-914
###################################################
# Panel 6(a)
betavalues <- betaCvMfit(PP)
crit.val <- qbeta((0.95)^(1 / 50), shape1 = betavalues[1],
    shape2 = betavalues[2])
par(mar=c(3,3,1,0.7), mgp=c(2,1,0))
hist(PP, breaks = 20, freq = FALSE, ylim = c(0, 50),
    col = "grey", main = "", xlab="Periodogram")
betafun <- function(x) dbeta(x, shape1 = betavalues[1],
    shape2 = betavalues[2])
curve(betafun, add = TRUE, lwd = 2)
abline(v = crit.val, lwd = 2)
par.mom <- betaCvMfit(PP, rob = FALSE, CvM = FALSE)
myf.mom <- function(x) dbeta(x, shape1 = par.mom[1],
    shape2 = par.mom[2])
curve(myf.mom, add = TRUE, lwd = 2, lty = 2)
crit.mom <- qbeta((0.95)^(1 / 50), shape1 = par.mom[1],
    shape2 = par.mom[2])
abline(v = crit.mom, lwd = 2, lty = 2)
par.rob <- betaCvMfit(PP, rob = TRUE, CvM = FALSE)
myf.rob <- function(x) dbeta(x, shape1 = par.rob[1],
    shape2 = par.rob[2])
curve(myf.rob, add = TRUE, lwd = 2, lty = 3)
crit.rob <- qbeta((0.95)^(1 / 50), shape1 = par.rob[1],
    shape2 = par.rob[2])
abline(v = crit.rob, lwd = 2, lty = 3)
legend("topright", lty = 1:3, legend = c("CvM", "Moments", "Robust moments"),
    bg = "white", lwd = 2)
box()


###################################################
### code chunk number 24: RobPer_Vignette.Rnw:918-923
###################################################
# Panel 6(b)
par(mar=c(3,3,1,0.1), mgp=c(2,1,0))
plot(1:50, PP, xlab = "Trial period", ylab = "Periodogram", main = "",
    type = "l")
abline(h = crit.val, lwd = 2)


