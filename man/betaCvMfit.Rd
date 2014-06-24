\encoding{utf8}
\name{betaCvMfit}
\alias{betaCvMfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Robust fit of a Beta distribution using CvM distance minimization
}
\description{
  Robustly fits a Beta distribution to \code{data} using \enc{Cramér}{Cramer}-von-Mises (CvM) distance minimization.
}
\usage{
betaCvMfit(data, CvM = TRUE, rob = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
     numeric vector: The sample, a Beta distribution is fitted to. 
}
  \item{CvM}{
     logical: If \code{FALSE} the \enc{Cramér}{Cramer}-von-Mises-distance is not minimized, but only moment estimates for the parameters of the Beta distribution are returned (see Details).
}
  \item{rob}{
     logical: If \code{TRUE}, mean and standard deviation are replaced by median and MAD when calculating moment estimates for the parameters of the Beta distribution (see Details). 
}
}
\details{
  \code{betaCvMfit} fits a Beta distribution to \code{data} by minimizing the \enc{Cramér}{Cramer}-von-Mises distance.
   Moment estimates of the parameters of the Beta distribution, clipped to positive values, are used as 
   starting values for the optimization process. They are calculated using
   \deqn{\hat a=-\frac{\bar x \cdot (-\bar x + \bar x^2 + \hat s^2)}{\hat s^2},}{a_est = - mean(x)*(-mean(x)+mean(x)^2+var(x))/var(x)}
   \deqn{\hat b= \frac{\hat a - \hat a \bar x}{\bar x}.}{b_est=(a_est-a_est*mean(x))/mean(x)}                                                                         
   
   These clipped moment estimates can be returned instead
   of CvM-fitted parameters setting \code{CvM = FALSE}.
   
   The \enc{Cramér}{Cramer}-von-Mises distance is defined as (see Clarke, McKinnon and Riley 2012)   
  \deqn{\frac 1n \sum_{i=1}^n \left(F(u_{(i)}) - 
  \frac{i-0.5}{n}\right)^2+ \frac{1}{12n^2},}{1/n*\sum(F(u[(i)])-(i-0.5)/n)^2+1/(12n^2)} where 
  \eqn{u_{(1)},\ldots,u_{(n)}}{u[(1)], \ldots, u[(n)]} is the ordered sample and \eqn{F} the 
  distribution function of Beta\eqn{(a,b)}{(a,b)}. 
}

\value{
numeric vector: Estimates for the Parameters \eqn{a,b}{a,b} of a Beta\eqn{(a,b)}{(a,b)} distribution with mean \eqn{a/(a+b)}{a/(a+b)}.
}
\author{
Anita M. Thieler, with contributions from Brenton R. Clarke.
}
\note{
Adapted from R-Code from Brenton R. Clarke to fit a Gamma distribution (see Clarke, McKinnon and Riley 2012) using \enc{Cramér}{Cramer}-von-Mises distance minimization. 
Used in Thieler et al. (2013). See also Thieler, Fried and Rathjens (2013).
}

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{RobPer-package}} for an example applying \code{betaCvMfit} to detect valid periods in a periodogram.
}

\references{
Clarke, B. R., McKinnon, P. L. and Riley, G. (2012): A Fast Robust Method for Fitting Gamma Distributions. Statistical Papers, 53 (4), 1001-1014

Thieler, A. M., Backes, M., Fried, R. and Rhode, W. (2013): Periodicity Detection in Irregularly Sampled Light Curves by Robust Regression and Outlier Detection. Statistical Analysis and Data Mining, 6 (1), 73-89

Thieler, A. M., Fried, R. and Rathjens, J. (2013): RobPer: An R Package to Calculate Periodograms for Light Curves Based on Robust Regression. Technical Report 2, Collaborative Research Center SFB 876, TU Dortmund University, available at \url{http://sfb876.tu-dortmund.de/PublicPublicationFiles/thieler_etal_2013b.pdf} 
}


\examples{
# data:
set.seed(12)
PP <- c(rbeta(45, shape1=4, shape2=15), runif(5, min=0.8, max=1))
hist(PP, freq=FALSE, breaks=30, ylim=c(0,7), xlab="Periodogram bar")

# true parameters:
myf.true <- function(x) dbeta(x, shape1=4, shape2=15)
curve(myf.true, add=TRUE, lwd=2)

# method of moments:
par.mom <- betaCvMfit(PP, rob=FALSE, CvM=FALSE)
myf.mom <- function(x) dbeta(x, shape1=par.mom[1], shape2=par.mom[2])
curve(myf.mom, add=TRUE, lwd=2, col="red")

# robust method of moments
par.rob <- betaCvMfit(PP, rob=TRUE, CvM=FALSE)
myf.rob <- function(x) dbeta(x, shape1=par.rob[1], shape2=par.rob[2])
curve(myf.rob, add=TRUE, lwd=2, col="blue")

# CvM distance minimization
par.CvM <- betaCvMfit(PP, rob=TRUE, CvM=TRUE)
myf.CvM <- function(x) dbeta(x, shape1=par.CvM[1], shape2=par.CvM[2])
curve(myf.CvM, add=TRUE, lwd=2, col="green")

# Searching for outliers...
abline(v=qbeta((0.95)^(1/50), shape1=par.CvM[1], shape2=par.CvM[2]), col="green")

legend("topright", fill=c("black", "green","blue", "red"),
    legend=c("true", "CvM", "robust moments", "moments"))
box()
}
