\name{lc_noise}
\alias{lc_noise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Noise and measurement accuracy generator for light curves
}
\description{
Generates measurement accuracies, a white noise component depending on them and a second (possibly power law, i.e. red) noise component which does not depend on the measurement accuracies.
For more details see \code{\link{tsgen}} or Thieler, Fried and Rathjens (2013). 
See \code{\link{RobPer-package}} for more information about light curves.
}
\usage{
lc_noise(tt, sig, SNR, redpart, alpha = 1.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tt}{
numeric vector: Observation times given.
}
  \item{sig}{
numeric vector of same length as \code{tt}: A given signal to which the noise will be added.
}
  \item{SNR}{
positive number: Defines the relation between signal and noise (see \code{\link{tsgen}} for Details).
}
  \item{redpart}{
numeric value in [0,1]: Proportion of the power law noise in noise components (see \code{\link{tsgen}} for Details).
}
  \item{alpha}{
numeric value: Power law index for the power law noise component (see \code{\link{tsgen}} for Details).
}
}


\value{
\item{y}{numeric vector: Observed values: signal + noise.}
\item{s}{numeric vector: Measurement accuracies related to the white noise component.}
}

\note{
A former version of this function is used in Thieler et al. (2013).
}

\references{
Thieler, A. M., Backes, M., Fried, R. and Rhode, W. (2013): Periodicity Detection in Irregularly Sampled Light Curves by Robust Regression and Outlier Detection. Statistical Analysis and Data Mining, 6 (1), 73-89

Thieler, A. M., Fried, R. and Rathjens, J. (2013): RobPer: An R Package to Calculate Periodograms for Light Curves Based on Robust Regression. Technical Report 2, Collaborative Research Center SFB 876, TU Dortmund University, available at \url{http://sfb876.tu-dortmund.de/PublicPublicationFiles/thieler_etal_2013b.pdf} 

}
\author{
Anita M. Thieler and Jonathan Rathjens
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Applied in \code{\link{tsgen}} (see there for an example), applies \code{\link{TK95_uneq}}.
}
