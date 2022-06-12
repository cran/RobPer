\name{tsgen}
\alias{tsgen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Artificial light curve generator
}
\description{
This function generates light curves (special time series) with unequally sampled observation times, different periodicities both in sampling and observed values, with white and power law (red) noise in the observed values and possibly disturbed observations.
See \code{\link{RobPer-package}} for more information about light curves and also Thieler, Fried and Rathjens (2016) for more details in general.
}
\usage{
tsgen(ttype, ytype, pf, redpart, s.outlier.fraction = 0, interval, npoints,
 ncycles, ps, SNR, alpha = 1.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ttype}{
character string: Specifying the sampling pattern. Possible choices are \code{"equi"} for equidistant sampling without gaps (unperiodic), \code{"unif"} for uniform non-equidistant unperiodic sampling, \code{"sine"} for sampling with periodic sine density, \code{"trian"} for sampling with periodic triangular density, both with period \eqn{p_s}{ps} (see Details and \code{\link{sampler}}). 
}
  \item{ytype}{
character string: Specifying the shape of the periodic fluctuation with period \eqn{p_f}{pf}. Possible choices are \code{"const"} for constantly being zero (so no periodicity), \code{"sine"} for a sine, \code{"trian"} for a periodic triangular function, \code{"peak"} for a peak function (see Details and \code{\link{signalgen}} for more details).
}
  \item{pf}{
positive number: Period \eqn{p_f}{pf} of the periodic fluctuation, argument of \code{signalgen} (see Details and \code{\link{signalgen}}).
}
  \item{redpart}{
numeric value in [0,1]: Proportion of the power law noise in noise components (see Details). The generated measurement accuracies \eqn{s_i}{s[i]} do not contain information about this noise component. 
}
  \item{s.outlier.fraction}{
numeric value in [0,1]: Fraction of measurement accuracies to be replaced by outliers. A value of 0 means that no measurement accuracy is replaced by an outlier (for more details see \code{\link{disturber}}). 
}
  \item{interval}{
logical: If \code{TRUE}, the observed values belonging to a random time interval of length 3\eqn{p_s}{ps} are replaced by atypical values (for more details see \code{\link{disturber}}).  
}
  \item{npoints}{
integer value: Defines the sample size \eqn{n}{n} for the generated light curve.
}
  \item{ncycles}{
integer value: number \eqn{n_s}{ns} of sampling cycles that is observed (see Details). 
}
  \item{ps}{
positive number: Sampling period \eqn{p_s}{ps}, influencing the sampling and how the light curve is disturbed (see Details and \code{\link{disturber}}).
}
  \item{SNR}{
positive number: Defines the relation between signal \eqn{y_f}{yf} and noise \eqn{y_w+y_r}{yw+yr} (see Details).
}
  \item{alpha}{
numeric value: Power law index \eqn{\alpha} for the power law noise component \eqn{y_r}{yr} (see Details). 
}
}
\details{
\code{tsgen} generates an artificial light curve consisting of observation times \eqn{t_1,\ldots,t_n}{t[1],\ldots,t[n]}, observation values \eqn{y_1,\ldots,y_n}{y[1],\ldots,y[n]} and measurement accuracies \eqn{s_1,\ldots,s_n}{s[1],\ldots,s[n]}. It calls several subfunctions (see there for details):

\code{\link{sampler}}  is used to sample observation times \eqn{t_1,\ldots,t_n}{t[1],\ldots,t[n]} in the interval \eqn{[0,n_s*p_s]}{0,ncycles*ps} with a possibly periodic sampling of period \eqn{p_s}{ps}.

\code{\link{signalgen}} generates periodically varying values \eqn{y_{f;1},\ldots,y_{f;n}}{yf[1],\ldots,yf[n]} at time points \eqn{t_1,\ldots,t_n}{t[1],\ldots,t[n]} with fluctuation period \eqn{p_f}{pf}.

\code{\link{lc_noise}}  samples measurement accuracies \eqn{s_1,\ldots,s_n}{s[1],\ldots,s[n]} 
from a Gamma(3,10)-distribution and a white noise component 
\eqn{y_{w;1},\ldots,y_{w;n}}{yw[1],\ldots,yw[n]} with from \eqn{\mathcal N(0,s_i^2)}{N(0,s[i]^2)} distributions. A second noise component \eqn{y_{r;1},\ldots,y_{r;n}}{yr[1],\ldots,yr[n]} does not depend on the \eqn{s_i}{s[i]}. It is generated as red noise, i.e. following a power law with power law index \eqn{\alpha}. For white noise choose \eqn{\alpha=0}, for flicker noise (pink noise) \eqn{\alpha=1}, for brown noise \eqn{\alpha=2}. The power law noise is generated using \code{\link{TK95_uneq}} and \code{\link{TK95}}. The noise components are scaled so that the variance of the \eqn{y_{r;i}}{yr[i]} has approximately the proportion \code{redpart} in the overall noise variance and that \code{SNR} is the ratio \eqn{var(y_f)/var(y_w+y_r)}{var(yf)/var(yw+yr)}. The observed values are set to \eqn{y_i= y_{f;i}+y_{w;i}+y_{r;i} \forall i}{y[i]=yf[i]+yw[i]+yr[i]}.

\code{\link{disturber}} disturbes the light curve replacing measurement accuracies 
\eqn{s_i}{s[i]} by outliers (if \code{s.outlier.fraction>0}) and observed values 
\eqn{y_i}{y[i]} by atypical values (if \code{interval=TRUE}). 
In case of \code{s.outlier.fraction=0} and \code{interval=FALSE}, the function returns all values unchanged.
}
\value{
  \item{tt}{numeric vector: Generated observation times\eqn{\ t_1,\ldots,t_n}.}    
  \item{y}{numeric vector: Generated observation values\eqn{\ y_1,\ldots,y_n}.} 
  \item{s}{numeric vector: Generated measurement accuracies\eqn{\ s_1,\ldots,s_n}.}                                       
}
\references{
Thieler, A. M., Backes, M., Fried, R. and Rhode, W. (2013): Periodicity Detection in Irregularly Sampled Light Curves by Robust Regression and Outlier Detection. Statistical Analysis and Data Mining, 6 (1), 73-89

Thieler, A. M., Fried, R. and Rathjens, J. (2016): RobPer: An R Package to Calculate Periodograms for Light Curves Based on Robust Regression. Journal of Statistical Software, 69 (9), 1-36, <doi:10.18637/jss.v069.i09>
}
\author{
Anita M. Thieler and Jonathan Rathjens
}
\note{
Note that the white noise components' variances are exactly \eqn{s_i^2}{s[i]^2}, so the \eqn{s_i}{s[i]} are no estimates, but true values. In this sense, the measurement accuracies of a generated light curve are more informative than for real light curves, where the measurement accuracies are estimates, see Thieler et al. (2013), where also a former version of this function is applied. 

To lower the informativity of the measurement accuracies, set \code{redpart} to a strictly positive value, possibly with \code{alpha=0} if no other noise components than white ones are required.
}


\seealso{
Applies \code{\link{sampler}}, \code{\link{signalgen}}, \code{\link{lc_noise}}, \code{\link{disturber}}, \code{\link{TK95}}, \code{\link{TK95_uneq}}.
}
\examples{
# Generate a light curve:
set.seed(22)
lightcurve<- tsgen(ttype="sine", ytype="peak" , pf=7, redpart=0.1, s.outlier.fraction=0, 
    interval=FALSE, npoints=200, ncycles=100, ps=5, SNR=3, alpha=0)

# Or do it step by step:
# First sampling observation times:
set.seed(22)
tt <- sampler(ttype="sine", npoints=200, ncycles=100, ps=5)

# show obviously irregular observation times as grey vertical bars on a red time line:
plot(tt, tt*0, type="n", axes=FALSE, xlab="Observation Times", ylab="")
    abline(v=tt, col="grey")
axis(1, pos=0, col="red", col.axis="red")

# Sampling period is 5, look at observation times modulo 10:
hist(tt\%\%5, xlab="Observation time modulo 5", 
    main="Sine Distribution for the phase (tt modulo 5)", freq=FALSE)
dsin <- function(tt) 0.2*(sin(2*pi*tt/5)+1)
curve(dsin, add=TRUE)

# Then generate periodic fluctuation
yf <- signalgen(tt, ytype="peak", pf=7)

plot(tt, yf, xlab="Observation Times", ylab="Periodic Fluctuation")
plot(tt\%\%7, yf, main="Phase Diagram (time modulo 7)", 
    xlab="Observation time modulo 7",  ylab="Periodic Fluctuation")

# Add noise and scale signal to the right SNR
temp <- lc_noise(tt,sig=yf, SNR=3, redpart=0.1, alpha=0)
y <- temp$y
s <- temp$s

# Plotting the light curve (vertical bars show measurement accuracies)
plot(tt, y, pch=16, cex=0.5, xlab="t", ylab="y", main="a Light Curve")
rect(tt, y+s, tt, y-s)

# The lightcurve has period 7:
plot(tt\%\%7, y, pch=16, cex=0.5, xlab="t", ylab="y", 
    main="Phase Diagram of a Light Curve")
rect(tt\%\%7, y+s, tt\%\%7, y-s)

# replace measurement accuracies by tiny outliers or include a peak
temp <- disturber(tt,y,s,ps=5, s.outlier.fraction=0, interval=FALSE)

# Phase diagram (observation times modulo 10)
plot(tt\%\%7, temp$y, pch=16, cex=0.5, xlab="t", ylab="y", 
    main="Phase Diagram of a Light Curve")
rect(tt\%\%7, temp$y+temp$s, tt\%\%7, temp$y-temp$s)

# The result is the same:
all(cbind(tt,temp$y,temp$s)==lightcurve)
}
