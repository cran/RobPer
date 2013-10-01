\encoding{utf8}
\name{TK95}
\alias{TK95}

\title{
Power law noise generator
}
\description{
Generates an equidistant time series of power law noise according to Timmer and \enc{König}{Koenig} (1995).
}
\usage{
TK95(N = 1000, alpha = 1.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
integer value: Length of the generated time series.
}
  \item{alpha}{
numeric value: Exponent of the power law. White noise has exponent 0, flicker noise (pink noise) has exponent 1, brown noise has exponent 2.
}
}

\value{
numeric vector: The generated time series.
}

\references{
Thieler, A. M., Backes, M., Fried, R. and Rhode, W. (2013): Periodicity Detection in Irregularly Sampled Light Curves by Robust Regression and Outlier Detection. Statistical Analysis and Data Mining, 6 (1), 73-89

Thieler, A. M., Fried, R. and Rathjens, J. (2013): RobPer: An R Package to Calculate Periodograms for Light Curves Based on Robust Regression. Technical Report 2, Collaborative Research Center SFB 876, TU Dortmund University, available at \url{http://sfb876.tu-dortmund.de/PublicPublicationFiles/thieler_etal_2013b.pdf} 

Timmer, J. and \enc{König}{Koenig}, M. (1995) On Generating Power Law Noise. Astronomy and Astrophysics, 300, 707-710
}



\author{
Anita M. Thieler with contributions of Uwe Ligges
}

\note{
This function is used in Thieler et al. (2013). See also Thieler, Fried and Rathjens (2013).
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Applied in \code{\link{tsgen}} by \code{\link{TK95_uneq}}.
}
\examples{
set.seed(31)
# Generate power law noise with exponent alpha=1.5:
y  <- TK95(N=2000, alpha=1.5)
tt <- seq(along=y)

# Show time series:
plot(tt,y, type="l", main="Power Law Noise", xlab="t", ylab="y")

# Plot Fourier periodogram with log-axes:
temp <- spectrum(y, plot=FALSE)
plot(log(temp$freq), log(temp$spec), main="log-log-Fourier periodogram",
     xlab="log(frequency)", ylab="log(periodogram)")

# A line with slope -alpha for comparison
abline(a=8, b=-1.5, col="red")
text(-2, 12, expression(alpha==1.5), col="red")
}
