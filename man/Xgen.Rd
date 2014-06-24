\name{Xgen}
\alias{Xgen}

\title{
Designmatrix generator
}
\description{
This function is used to create the designmatrices needed in \code{RobPer} to fit periodic functions. See \code{RobPer} or Thieler, Fried and Rathjens (2013) for Details.
}
\usage{
Xgen(tt, n, s, pp, design, steps = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tt}{
real-valued vector of length \code{n}: Observation times. 
}
  \item{n}{
integer: Sample size.
}
  \item{s}{
numeric vector of length \code{n}: Measurement errors to perform weighted regression, else set \code{s=rep(1,n)}.
}
  \item{pp}{
positive number: Trial period.
}
  \item{design}{
character string: Shape of the periodic function to be fitted, possible choices are \code{"step"}, \code{"sine"}, \code{"fourier(2)"}, \code{"fourier(3)"}, \code{"splines"} (see \code{RobPer} for details) and  \code{"stepB"}. The latter generates a step function like \code{"step"}, but with opposite jumping time points. This is needed for RobPer with \code{model="2step"} (see \code{RobPer}).
}
  \item{steps}{
Number of steps for the step functions 
}
}

\value{
numeric matrix: Designmatrix.
}
\references{
Thieler, A. M., Backes, M., Fried, R. and Rhode, W. (2013): Periodicity Detection in Irregularly Sampled Light Curves by Robust Regression and Outlier Detection. Statistical Analysis and Data Mining, 6 (1), 73-89

Thieler, A. M., Fried, R. and Rathjens, J. (2013): RobPer: An R Package to Calculate Periodograms for Light Curves Based on Robust Regression. Technical Report 2, Collaborative Research Center SFB 876, TU Dortmund University, available at \url{http://sfb876.tu-dortmund.de/PublicPublicationFiles/thieler_etal_2013b.pdf} 

}
\author{
Anita M. Thieler and Jonathan Rathjens
}
\note{
A former version of this function is used in Thieler et al. (2013).
}


\seealso{
Applied in \code{\link{RobPer}}, see \code{\link{FastTau}} for an example.
}
