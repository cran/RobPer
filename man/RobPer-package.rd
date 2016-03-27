\encoding{latin1}
\name{RobPer-package}
\alias{RobPer-package}
\docType{package}
\title{
The RobPer-package
}
\description{
Calculates periodograms based on (robustly) fitting periodic functions to light curves and other irregulary observed time series and detects high periodogram bars.
}
\details{
\tabular{ll}{
Package: \tab RobPer\cr
Type: \tab Package\cr
Version: \tab 1.2.2\cr
Date: \tab 2016-03-27\cr
License: \tab GPL-3 \cr
}
Light curves occur in astroparticle physics and are irregularly sampled times series
\eqn{(t_i, y_i)_{i=1,\ldots,n}}{(t[i], y[i]), i=1,\ldots,n,}  or \eqn{(t_i, y_i, s_i)_{i=1,\ldots,n}}{(t[i], y[i], s[i]), i=1,\ldots,n,} consisting of unequally 
spaced observation times \eqn{t_1, \ldots, t_n}{t[1],\ldots,t[n]}, observed values
\eqn{y_1, \ldots, y_n}{y[1],\ldots,y[n]} and possibly measurement accuracies \eqn{s_1, \ldots, s_n}{s[1],\ldots,s[n]}.
The pattern of the observation times \eqn{t_i}{t[i]} may be periodic with sampling period \eqn{p_s}{ps}. 
The observed values \eqn{y_i}{y[i]} may possibly contain a periodic fluctuation \eqn{y_{f;i}}{yf[i]} with 
fluctuation period \eqn{p_f}{pf}.  One is interested in finding \eqn{p_f}{pf}. The measurement 
accuracies \eqn{s_i}{s[i]} give information about how precise the \eqn{y_i}{y[i]} were measured. 
They can be interpreted as estimates for the standard deviations of the observed values.
For more details see Thieler et al. (2013) or Thieler, Fried and Rathjens (2016). 

This package includes three main functions: \code{RobPer} calculates the periodogram, possibly taking into account measurement accuracies. With \code{betaCvMfit}, outlying periodogram bars (indicating a period) can be detected. This function bases on robustly fitting a distribution using \enc{Cramér}{Cramer}-von-Mises (CvM) distance minimization (see also Clarke, McKinnon and Riley 2012).  The function \code{tsgen} can be used to generate artificial light curves. For more details about the implementation see Thieler, Fried and Rathjens (2016).

A preliminary version of this package is used in Thieler et al. (2013). The \code{FastS}-function  
and the \code{FastTau}-function
presented here are slightly changed versions of R-Code published in Salibian-Barrera and Yohai (2006)
and Salibian-Barrera, Willems and Zamar (2008).

The financial support of the DFG (SFB 876 "Providing Information by Resource-Constrained Data Analysis", project C3, and GK 1032 "Statistische Modellbildung") is gratefully acknowledged. We thank the ITMC at TU Dortmund University for providing computer resources on LiDO.
}
\author{
Anita M. Thieler, Jonathan Rathjens and Roland Fried, with contributions from Brenton R. Clarke (see  \code{\link{betaCvMfit}}),  Matias Salibian-Barrera, Gert Willems and Victor Yohai (see \code{\link{FastS}} and \code{\link{FastTau}}) and Uwe Ligges (see \code{\link{TK95}}).

Maintainer: Jonathan Rathjens <jonathan.rathjens@tu-dortmund.de>
}
\references{
Clarke, B. R., McKinnon, P. L. and Riley, G. (2012): A Fast Robust Method for Fitting Gamma Distributions. Statistical Papers, 53 (4), 1001-1014

Salibian-Barrera, M. and Yohai, V. (2006): A Fast Algorithm for 
S-Regression Estimates. Journal of Computational and Graphical 
Statistics, 15 (2), 414-427 

Salibian-Barrera, M., Willems, G. and Zamar, R. (2008): The Fast-tau Estimator for Regression. Journal of Computational and Graphical Statistics, 17 (3), 659-682

Thieler, A. M., Backes, M., Fried, R. and Rhode, W. (2013): Periodicity Detection in Irregularly Sampled Light Curves by Robust Regression and Outlier Detection. Statistical Analysis and Data Mining, 6 (1), 73-89

Thieler, A. M., Fried, R. and Rathjens, J. (2016): RobPer: An R Package to Calculate Periodograms for Light Curves Based on Robust Regression. Journal of Statistical Software, 69 (9), 1-36, <doi:10.18637/jss.v069.i09>
}

\examples{
# Generate a disturbed light curve:
set.seed(22)
lightcurve <- tsgen(ttype="sine",ytype="peak" , pf=7, redpart=0.1, s.outlier.fraction=0.1,
    interval=TRUE, npoints=200, ncycles=25, ps=20, SNR=3, alpha=0)

# Plotting the light curve (vertical bars show measurement accuracies)
plot(lightcurve[,1], lightcurve[,2], pch=16, cex=0.5, xlab="t", ylab="y", 
    main="a Light Curve")
rect(lightcurve[,1], lightcurve[,2]+lightcurve[,3], lightcurve[,1], 
    lightcurve[,2]-lightcurve[,3])

# The lightcurve has a period of 7:
plot(lightcurve[,1]\%\%7, lightcurve[,2], pch=16, cex=0.5, xlab="t", ylab="y",
    main="Phase Diagram of a Light Curve")
rect(lightcurve[,1]\%\%7, lightcurve[,2]+lightcurve[,3], lightcurve[,1]\%\%7, 
    lightcurve[,2]-lightcurve[,3])

# Calculate a periodogram of a light curve:
PP <- RobPer(lightcurve, model="splines", regression="huber", weighting=FALSE, 
    var1=FALSE, periods=1:50)

# Searching for extremely high periodogram bars:
betavalues <- betaCvMfit(PP)
crit.val <- qbeta((0.95)^(1/50),shape1=betavalues[1], shape2=betavalues[2])

hist(PP, breaks=20, freq=FALSE, ylim=c(0,100), xlim=c(0,0.08), col=8, main ="")
betafun <- function(x) dbeta(x, shape1=betavalues[1], shape2=betavalues[2])
curve(betafun, add=TRUE, lwd=2)
abline(v=crit.val, lwd=2)

# alternatives for fitting beta distributions:
# method of moments:
par.mom <- betaCvMfit(PP, rob=FALSE, CvM=FALSE)
myf.mom <- function(x) dbeta(x, shape1=par.mom[1], shape2=par.mom[2])
curve(myf.mom, add=TRUE, lwd=2, col="red")
crit.mom <- qbeta((0.95)^(1/50),shape1=par.mom[1], shape2=par.mom[2])
abline(v=crit.mom, lwd=2, col="red")

# robust method of moments
par.rob <- betaCvMfit(PP, rob=TRUE, CvM=FALSE)
myf.rob <- function(x) dbeta(x, shape1=par.rob[1], shape2=par.rob[2])
curve(myf.rob, add=TRUE, lwd=2, col="blue")
crit.rob <- qbeta((0.95)^(1/50),shape1=par.rob[1], shape2=par.rob[2])
abline(v=crit.rob, lwd=2, col="blue")

legend("topright", fill=c("black","red","blue"), 
    legend=c("CvM", "moments", "robust moments"), bg="white")
box()

# Detect fluctuation period:
plot(1:50, PP, xlab="Trial Period", ylab="Periodogram", type="l", 
    main="Periodogram fitting periodic splines using M-regression (Huber function)")
abline(h=crit.val, lwd=2)
text(c(7,14), PP[c(7,14)], c(7,14), adj=1, pos=4)
axis(1, at=7, labels=expression(p[f]==7))

# Comparison with non-robust periodogram
# (see package vignette, section 5.1 for further graphical analysis)
PP2 <- RobPer(lightcurve, model="splines", regression="L2",
    weighting=FALSE, var1=FALSE, periods=1:50)
betavalues2 <- betaCvMfit(PP2)
crit.val2   <- qbeta((0.95)^(1/50),shape1=betavalues2[1], shape2=betavalues2[2])

plot(1:50, PP2, xlab="Trial Period", ylab="Periodogram", type="l", 
    main="Periodogram fitting periodic splines using L2-regression")
abline(h=crit.val2, lwd=2)
}
