\name{disturber}
\alias{disturber}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Disturbing light curve data
}
\description{
Disturbes a light curve replacing measurement accuracies by outliers and/or observed values by atypical values.
See \code{\link{RobPer-package}} for more information about light curves.
}
\usage{
disturber(tt, y, s, ps, s.outlier.fraction = 0, interval)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tt}{
numeric vector: Observation times \eqn{t_1,\ldots,t_n}{} (see Details).
}
  \item{y}{
numeric vector: Observed values \eqn{y_1,\ldots,y_n}{} (see Details).
}
  \item{s}{
numeric vector: Measurement accuracies \eqn{s_1,\ldots,s_n}{} (see Details).
}
  \item{ps}{
positive value: Sampling period\eqn{\ p_s.}{.} Indirectly defines the length of the time interval, in which observed values \eqn{y_i}{} are replaced by atypical values (see Details).
}
  \item{s.outlier.fraction}{
numeric value in [0,1]: Defines the proportion of measurement accuracies that is replaced by outliers (see Details). A value of 0 means that no measurement accuracy is replaced by an outlier.
}
  \item{interval}{
logical: If \code{TRUE}, the observed values belonging to a random time interval of length 3\eqn{p_s}{ps} are replaced by atypical values (see Details). If \code{TRUE} and the light curve is shorter than \eqn{3p_s}{ps}, the function will stop with an error message. 
}
}
\details{
This function disturbes the light curve \eqn{(t_i,y_i,s_i)_{i=1,\ldots,n}}{(t[i], y[i], s[i]), i=1,\ldots,n,} given. It randomly chooses a proportion of \code{s.outlier.fraction} measurement accuracies \eqn{s_i}{s} and replaces them by \eqn{0.5\min(s_1,\ldots,s_n)}{0.5*min(s)}. In case of \code{interval=TRUE} a time interval \eqn{[t_{start},t_{start}+3p_s]}{[t.start, t.start+3ps]} within the intervall 
\eqn{[t_1,t_n]}{[t[1], t[n]]} is randomly chosen and all observed values belonging to this time interval are replaced by a peak function:
\deqn{y_i^{changed} = 6 \ \tilde y_{0.9}\  \frac{d_{\mathcal N(t_{start}+1.5p_s, p_s^2)}(t_i) }{ d_{\mathcal N(0,p_s^2)}(0)} \quad \forall \ i \ : \ t_i\in[t_{start}, t_{start}+3p_s],}{y[i]=6*quantile(y, 0.9)*dnorm(t[i], mean=t.start+1.5ps, sd=ps)/dnorm(0, sd=ps), for i : t[i] in [t.start, t.start+3ps],} 
where \eqn{d_{\mathcal N(a,b^2)}(x)}{dnorm(x, mean=a, sd=b)} denotes the density of a normal distribution with mean \eqn{a} and variance \eqn{b^2} at \eqn{x}.

In case of \code{s.outlier.fraction=0} and \code{interval=FALSE}, \code{y} and \code{s} are returned unchanged.
}
\value{
  \item{y}{numeric vector: New \eqn{y_i}{y}-values, partly different from the old ones if \code{interval=TRUE} (see Details).}
  \item{s}{numeric vector: New \eqn{s_i}{s}-values, partly different from the old ones if \code{s.outlier.fraction>0} (see Details).}
}
\references{
Thieler, A. M., Backes, M., Fried, R. and Rhode, W. (2013): Periodicity Detection in Irregularly Sampled Light Curves by Robust Regression and Outlier Detection. Statistical Analysis and Data Mining, 6 (1), 73-89

Thieler, A. M., Fried, R. and Rathjens, J. (2013): RobPer: An R Package to Calculate Periodograms for Light Curves Based on Robust Regression. Technical Report 2, Collaborative Research Center SFB 876, TU Dortmund University, available at \url{http://sfb876.tu-dortmund.de/PublicPublicationFiles/thieler_etal_2013b.pdf} 

}
\author{
Anita M. Thieler
}
\note{
A former version of this function is used in Thieler et al. (2013). See also Thieler, Fried and Rathjens (2013).
}


\seealso{
Applied in \code{\link{tsgen}} (see there for example).
}


