\encoding{utf8}
\name{TK95_uneq}
\alias{TK95_uneq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power law noise generator for unequally sampled observation times
}
\description{
Generates power law noise using \code{\link{TK95}} according to Timmer and \enc{König}{Koenig} (1995), with modifications proposed in Uttley, McHardy and Papadakis (2002) for given irregular observation times.

}
\usage{
TK95_uneq(tt, alpha = 1.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tt}{
numeric vector: Observation times given.
}
  \item{alpha}{
numeric value: exponent of the power law. White noise has exponent 0, flicker noise (pink noise) has exponent 1, brown noise has exponent 2.
}
}


\value{
numeric vector: Noise values related to the observation times.
}

\note{
This function is applied in Thieler et al. (2013). See also Thieler, Fried and Rathjens (2013).
}
\references{
Thieler, A. M., Backes, M., Fried, R. and Rhode, W. (2013): Periodicity Detection in Irregularly Sampled Light Curves by Robust Regression and Outlier Detection. Statistical Analysis and Data Mining, 6 (1), 73-89

Thieler, A. M., Fried, R. and Rathjens, J. (2013): RobPer: An R Package to Calculate Periodograms for Light Curves Based on Robust Regression. Technical Report 2, Collaborative Research Center SFB 876, TU Dortmund University, available at \url{http://sfb876.tu-dortmund.de/PublicPublicationFiles/thieler_etal_2013b.pdf} 

Timmer, J. and \enc{König}{Koenig}, M. (1995) On Generating Power Law Noise. Astronomy and Astrophysics, 300, 707-710

Uttley, P., McHardy, I. M. and Papadakis, I. E. (2002) Measuring the Broad-Band Power Spectra of Active Galactic Nuclei with RXTE. Monthly Notices of the Royal Astronomical Society, 332 (1), 231-250
}
\author{
Anita M. Thieler
}


\seealso{
Applies \code{\link{TK95}}, applied in \code{\link{tsgen}}.
}

\examples{
# Compare with example in TK95 to see that the power law is much more clear in
# equally sampled data!
set.seed(31)
# Generate power law noise with exponent alpha=1.5:
tt <- sampler(ttype="unif", ps=1, ncycles=2000, npoints=2000)
y <- TK95_uneq(tt, alpha=1.5)

# Show time series:
plot(tt,y, type="l", main="Irregular Power Law Noise", xlab="t", ylab="y")

# Plot Lomb-Scargle periodogram with log-axes:
temp <- RobPer(cbind(tt,y,1), weighting=FALSE, model="sine", regression="L2",
    periods=2000/seq(2, 1000, 2))
plot(log(seq(2, 1000, 2)/2000), log(temp), main="log-log-Fourier periodogram",
    xlab="log(frequency)", ylab="log(periodogram)")
title(main= "Power Law not so obvious", cex.main=0.8, line=0.5)

# A line with slope -alpha for comparison
abline(a=-10, b=-1.5, col="red")
text(-5, -1.5, expression(alpha==1.5), col="red")
}
