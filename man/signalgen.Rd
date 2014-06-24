\name{signalgen}
\alias{signalgen}

\title{
Generator for periodic signal in a light curve
}
\description{
Calculates periodically varying values for given observation times.

}
\usage{
signalgen(tt, ytype, pf = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tt}{
numeric vector: Observation times \eqn{t_1,\ldots,t_n}{} (see Details).
}
  \item{ytype}{
character string: Specifying the shape of the periodic fluctuation (see Details). Possible choices are \code{"const"}, \code{"sine"}, \code{"trian"},\code{"peak"}.

}
  \item{pf}{
positive numeric value: Fluctuation period\eqn{ p_f.}{.}
}
}
\details{
The values \eqn{y_{f;1},\ldots,y_{f;n}}{yf[1],\ldots,yf[n]} with fluctuation period \eqn{p_f}{pf} and related to observation times \eqn{t_1,\ldots,t_n}{t[1],\ldots,t[n]} are generated using 
\deqn{y_{f;i}=f\left(\frac{t_i}{p_f}\right), i=1,\ldots,n. }{yf[i]=f(t[i]/pf), i=1,\ldots,n.}
Depending on \code{ytype} (see above), \eqn{f} is defined as:

\deqn{ f_{const}(t) = 0,}{f.const(t)=0,}

\deqn{f_{sine}(t)= \sin\left(\frac{2\pi t}{p_f}\right),}{f.sine(t)=sin(2\pi t/pf),}

\deqn{f_{trian}(t)= 3\varphi_{1}(t), \quad 0\leq \varphi_{1}(t)\leq\frac{2}{3},}{f.trian(t)=3*phase(t) if 0 \le phase(t) \le 2/3,}

\deqn{f_{trian}(t)= 6-6\varphi_{1}(t),\quad \frac{2}{3}<\varphi_{1}(t)\leq 1,}{f.trian(t)=6-6*phase(t) if 2/3 < phase(t) \le 1,}

\deqn{f_{peak}(t)= 9\exp\left(-3p_f^2\left(\varphi_{1}(t)-\frac 23\right)^2\right),\quad  0\leq \varphi_{1}(t)\leq\frac{2}{3},}{f.peak(t)=9exp(-3pf^2*(phase(t)-2/3)^2) if 0 \le phase(t) \le 2/3,}

\deqn{f_{peak}(t)= 9\exp\left(-12p_f^2\left(\varphi_{1}(t)-\frac 23\right)^2\right),\quad \frac{2}{3}<\varphi_{1}(t)\leq 1,}{f.peak(t)=9exp(-12pf^2*(phase(t)-2/3)^2) if 2/3 < phase(t) \le 1,}

with \eqn{\varphi_1(t) = t mod1 = (t-\lfloor t/p_f \rfloor p_f)/p_f}{phase(t) = t mod 1} = \code{(t\%\%1)/pf}. \eqn{f_{const}}{f.const} means that there is no (periodic) fluctuation, \eqn{f_{sine}}{f.sine} defines a sine function, \eqn{f_{trian}}{f.trian} defines a triangular shaped periodic function and \eqn{f_{peak}}{f.peak} a periodically repeating peak.
}
\value{
numeric vector: Values \eqn{y_{f;1},\ldots,y_{f;n}}{yf[1],\ldots,yf[n]}.
}

\note{This function is used in Thieler et al. (2013). See also Thieler, Fried and Rathjens (2013).}

\references{
Thieler, A. M., Backes, M., Fried, R. and Rhode, W. (2013): Periodicity Detection in Irregularly Sampled Light Curves by Robust Regression and Outlier Detection. Statistical Analysis and Data Mining, 6 (1), 73-89

Thieler, A. M., Fried, R. and Rathjens, J. (2013): RobPer: An R Package to Calculate Periodograms for Light Curves Based on Robust Regression. Technical Report 2, Collaborative Research Center SFB 876, TU Dortmund University, available at \url{http://sfb876.tu-dortmund.de/PublicPublicationFiles/thieler_etal_2013b.pdf} 
}
\author{
Anita M. Thieler and Jonathan Rathjens
}

\seealso{
Applied in \code{\link{tsgen}} (see there for an example).
}
