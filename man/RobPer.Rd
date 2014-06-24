\encoding{utf8}
\name{RobPer}
\alias{RobPer}

\title{
Periodogram based on (robustly) fitting a periodic function to a light curve
}
\description{
Calculates a periodogram by fitting a periodic function to a light curve, using a possibly robust regression technique and possibly taking into account measurement accuracies.  
See \code{\link{RobPer-package}} for more information about light curves. For a lot of more details see Thieler, Fried and Rathjens (2013) and Thieler et al. (2013).
}
\usage{
RobPer(ts, weighting, periods, regression, model, steps = 10, tol = 1e-03,
 var1 = weighting, genoudcontrol = list(pop.size = 50, max.generations = 50,
 wait.generations = 5), LTSopt =TRUE, 
 taucontrol = list(N = 100, kk = 2, tt = 5, rr = 2, approximate = FALSE),
 Scontrol=list(N = ifelse(weighting,200,50), kk = 2, tt = 5, b=.5, cc = 1.547,
 seed = NULL) )
}

\arguments{
  \item{ts}{
dataframe or matrix with three (or two) numeric columns containing the light curve to be analyzed: 
observation times (first column), observed values (second column), measurement accuracies (thirs column). 
If it is intended to calculate the periodogram of a time series without measurement accuracies (\code{weighting=FALSE}), 
the third column may be omitted.
}
  \item{weighting}{
logical: Should measurement accuracies be taken into account performing weighted regression?
}
  \item{periods}{
vector of positive numeric values: Trial periods.
}
  \item{regression}{
character string specifying the regression method used: Possible choices are
\code{"L2"} (least squares regression using the R-function \code{lm}, package \code{stats}), 
\code{"L1"} (least absolute deviation regression, using the R-function \code{rq}, package \code{quantreg}),
\code{"LTS"} (least trimmed squares regression, using the R-function \code{ltsReg}, package \code{robustbase}),
\code{"huber"} (M-regression using the Huber function),
\code{"bisquare"} (M-regression using the bisquare function),
\code{"S"} (S-regression using adapted code from Salibian-Barrera and Yohai 2006, see \code{\link{FastS}}),
\code{"tau"} (tau-regression using adapted code from Salibian-Barrera, Willems and Zamar 2008, see \code{\link{FastTau}}).
}
  \item{model}{
character string specifying the periodic function fitted to the light curve: Possible choices are
\code{"step"} (periodic step function),
\code{"2step"} (two overlapping periodic step functions, see Details),
\code{"sine"}  (sine function),
\code{"fourier(2)"} and \code{"fourier(3)"} (Fourier series of second or third degree),
\code{"splines"} (periodic spline function with four B-splines per cycle, generated using \code{spline.des}, package \code{splines}).
}
  \item{steps}{
integer value: Number of steps per cycle for the periodic step function(s).
}
  \item{tol}{
(small) positive number: Precision for convergence criteria. Used in case of \code{regression="huber"} or \code{"bisquare"} or if \code{regression="LTS"} and \code{LTSopt=TRUE}.
}
  \item{var1}{
logical: Should variance estimate be set to 1 in case of weighted M-regression?
}
  \item{genoudcontrol}{
list of three integers \code{pop.size}, \code{max.generations}, \code{wait.generations}: Control parameters for the R-function \code{genoud}, package \code{rgenoud}, see Details and Mebane Jr. and Sekhon (2011). Used in case of \code{regression="bisquare"} or if \code{regression="LTS"} and \code{LTSopt=TRUE}.
}
  \item{LTSopt}{
logical: In case of LTS-regression, should regression result of \code{ltsReg} be optimized using the R-function \code{genoud}, package \code{rgenoud}?
}
  \item{taucontrol}{
list of four integer values \code{N}, \code{kk}, \code{tt}, \code{rr} and one logical \code{approximate}: Control parameters for the R-function \code{FastTau}. For more details see \code{\link{FastTau}} and Salibian-Barrera, Willems and Zamar (2008).
}
  \item{Scontrol}{
list of three integers \code{N}, \code{kk} and \code{tt}, two positive numbers \code{b} and \code{cc} and another integer \code{seed}: Control parameters for the R-function \code{FastS}. For more details see \code{\link{FastS}} and Salibian-Barrera and Yohai (2006). Please notice that the further Scontrol entry \code{int} expected by \code{FastS} is automatically set to \code{FALSE} in order to let \code{RobPer} work properly.
} 
}
\details{
For each trial period, a periodic function (defined by \code{model}) is fitted to the light curve using regression technique \code{regression}. The periodogram bar is the coefficient of determination. In case of \code{model="2step"}, two different step functions with opposed jumping times are fitted separately and the periodogram bar is the mean of both coefficients of determination. For a lot of more details see Thieler, Fried and Rathjens (2013) and Thieler et al. (2013).
}
\value{
numeric vector: Periodogram bars related to the trial periods.
}
\references{
Leahy, D. A., Darbro, W., Elsner, R. F., Weisskopf, M. C., Kahn, S., Sutherland, P. G. and Grindlay, J. E. (1983): On Searches for Pulsed Emission with Application to Four Globular Cluster X-ray Sources-NGC 1851, 6441, 6624, and 6712. The Astrophysical Journal, 266 (1), 160-170

Mebane Jr., W. R. and Sekhon, J. S. (2011): Genetic Optimization Using Derivatives: The rgenoud Package for R. Journal of Statistical Software, 42 (11), 1-26

Salibian-Barrera, M. and Yohai, V. (2006): A Fast Algorithm for S-Regression Estimates. Journal of Computational and Graphical Statistics, 15 (2), 414-427 

Salibian-Barrera, M., Willems, G. and Zamar, R. (2008): The Fast-tau Estimator for Regression. Journal of Computational and Graphical Statistics, 17 (3), 659-682

Stellingwerf, R. F. (1978): Period Determination Using Phase Dispersion Minimization. The Astrophysical Journal, 224, 953-960

Thieler, A. M., Backes, M., Fried, R. and Rhode, W. (2013): Periodicity Detection in Irregularly Sampled Light Curves by Robust Regression and Outlier Detection. Statistical Analysis and Data Mining, 6 (1), 73-89

Thieler, A. M., Fried, R. and Rathjens, J. (2013): RobPer: An R Package to Calculate Periodograms for Light Curves Based on Robust Regression. Technical Report 2, Collaborative Research Center SFB 876, TU Dortmund University, available at \url{http://sfb876.tu-dortmund.de/PublicPublicationFiles/thieler_etal_2013b.pdf} 

Zechmeister, M. and \enc{Kürster}{Kuerster}, M. (2009): The Generalised Lomb-Scargle Periodogram. A New Formalism for the Floating-Mean and Keplerian Periodograms. Astronomy and Astrophysics, 496 (2), 577-584
}
\author{
Anita M. Thieler, Jonathan Rathjens and Roland Fried
}
\note{
Performing \code{weighting = FALSE}, \code{regression="L2"}, \code{model="sine"} on a equidistantly sampled time series is equivalent to calculating the standard periodogram of Fourier analysis, see Example.

Performing \code{regression="L2"}, \code{model="sine"} is equivalent to calculating a Generalized Lomb-Scargle periodogram (see Zechmeister and \enc{Kürster}{Kuerster} 2009).

Performing \code{regression="L2"}, \code{model="step"} is equivalent to calculating an Epoch Folding (Leahy et al. 1983) or Anaysis of Variance (Schwarzenberg-Czerny 1989) periodogram.

Performing \code{regression="L2"}, \code{model="2step"} is equivalent to calculating a Phase Dispersion Minimization periodogram (Stellingwerf 1978).

A former version of this function is used in Thieler et al. (2013). For more equivalences see there. 
}

\examples{
# For more examples see RobPer-package and TK95_uneq!

# Example to show the equivalence between the periodogram from Fourier analysis
# and the Lomb-Scargle periodogram in case of equidistant sampling and equal weighting:
set.seed(7)
n <- 120
# equidistant time series:
zr <- tsgen(ttype="equi", ytype="const", pf=1, redpart= 0, s.outlier.fraction=0.2, 
    interval=FALSE, npoints=n, ncycles=n, ps=1, SNR=1, alpha=1.5)
# periodogram of Fourier analysis
PP_konv <- spec.pgram(zr[,2], taper = 0,pad = 0, fast = FALSE, demean = TRUE,
    detrend = TRUE, plot = TRUE)
# Lomb-Scargle periodogram - Note: Due to the regression ansatz,
# RobPer is not able to compute period 2 in this case.
PP_new <- RobPer(ts=zr, weighting=FALSE, periods=1/PP_konv$freq,
    regression="L2", model="sine")
plot(PP_konv$freq, PP_konv$spec, ylab="periodogram", xlab="frequency",
    main="Comparison of RobPer(...regression='LS', model='sine') and spec.pgram")
points(PP_konv$freq, PP_new*var(zr[,2])*n/2, type="l")
legend("top",lty=c(1,0), pch=c(-5,1), legend=c("RobPer*var(y)*n/2", "spec.pgram"))
# Due to different ways of computation, the scaled periodograms are not exactly
# identical, but show very similar behavior.
}

\seealso{
Applies \code{\link{FastS}} and \code{\link{FastTau}}, \code{\link{Xgen}}, examples in \code{\link{RobPer-package}} and \code{\link{TK95_uneq}}. 
}
