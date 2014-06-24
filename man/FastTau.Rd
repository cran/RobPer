\name{FastTau}
\alias{FastTau}

\title{
Tau-Regression using the Fast-tau-Algorithm
}
\description{
Performs tau-Regression using the Fast-tau-Algorithm.
}
\usage{
FastTau(x, y, taucontrol = list(N = 500, kk = 2, tt = 5, rr = 2, approximate = 0),
 beta_gamma)
}

\arguments{
  \item{x}{
numeric \eqn{(n\times p)}{(n,p)}-matrix: Designmatrix.
}
  \item{y}{
numeric vector: \eqn{n} observations.
}
  \item{taucontrol}{
list of four integer and one logical value: Control parameters (see Details).
}
  \item{beta_gamma}{
numeric vector: Specifies one parameter candidate of length \eqn{p} (see Details).
}
}
\details{
The Fast-tau-Algorithm to efficiently perform tau-Regression was published by Salibian-Barrera, Willems and Zamar (2008).
It bases on starting with a set of \code{N} parameter candidates,
locally optimizing them using \code{kk} iterations, 
then optimizing the \code{tt} best candidates to convergence and finally choosing the best parameter candidate. Since calculation of the objective value is computationally expensive, it is possible to approximate it with \code{rr} iteration steps when choosing \code{approximate=TRUE}. For more details see Salibian-Barrera, Willems and Zamar (2008).

The R-function \code{FastTau} used in \code{RobPer} is a slightly changed version of the R-code published in Salibian-Barrera, Willems and Zamar (2008). It was changed in order to work more efficiently, especially when fitting step functions, and to specify one parameter candidate in advance. For details see Thieler, Fried and Rathjens (2013).
}
\value{
\item{beta}{numeric vector: Fitted parameter vector.}
\item{scale}{numeric: Value of the objective function}
}

\references{
Salibian-Barrera, M., Willems, G. and Zamar, R. (2008): The Fast-tau Estimator for Regression. Journal of Computational and Graphical Statistics, 17 (3), 659-682

Thieler, A. M., Fried, R. and Rathjens, J. (2013): RobPer: An R Package to Calculate Periodograms for Light Curves Based on Robust Regression. Technical Report 2, Collaborative Research Center SFB 876, TU Dortmund University, available at \url{http://sfb876.tu-dortmund.de/PublicPublicationFiles/thieler_etal_2013b.pdf} 

}

\author{
Matias Salibian-Barrera and Gert Willems, modified by Anita M. Thieler
}


\seealso{
Applied in \code{\link{RobPer}}.
}
\examples{
set.seed(22)
# Generate a disturbed light curve
lightcurve <- tsgen(ttype="unif",ytype="sine" , pf=7, redpart=0.1, interval=TRUE,
    npoints=100, ncycles=10, ps=7, SNR=4, alpha=0)
tt <- lightcurve[,1]
y  <- lightcurve[,2]
s  <- rep(1,100)  # unweighted regression

plot(tt, y, type="l", main="Fitting a sine to a disturbed lightcurve")

# Fit the true model (a sine of period 7)... designmatrix:
X <- Xgen(tt, n=100, s, pp=7, design="sine")
# Robust tau-fit:
beta_FastTau <- FastTau(X, y)$beta
# Robust S-fit:
beta_FastS <- FastS(X, y)$beta
# Least squares fit:
beta_lm <- lm(y~0+X)$coeff

# Plot:
sin7_fun <- function(t, beta) beta[1]+ beta[2]*sin(t*2*pi/7)+ beta[3]*cos(t*2*pi/7)
sin_FastTau <- function(t) sin7_fun(t, beta_FastTau)
sin_FastS <- function(t) sin7_fun(t, beta_FastS)
sin_lm <- function(t) sin7_fun(t, beta_lm)
curve(sin_FastTau, col="green", add=TRUE)
curve(sin_FastS, col="blue", add=TRUE, lty=2)
curve(sin_lm, col="red", add=TRUE)

legend("topleft", fill=c("black", "red", "green", "blue"),
    legend=c("Light Curve (disturbed)", "Least Squares Fit", "FastTau Fit", "FastS Fit"))
}
