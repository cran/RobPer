\name{sampler}
\alias{sampler}

\title{
Generator for irregularly sampled observation times 
}
\description{
Generates irregularly sampled observation times with a periodic sampling pattern
}
\usage{
sampler(ttype, npoints, ncycles, ps = 1)
}
\arguments{
  \item{ttype}{
character string: Specifying the sampling pattern. Possible options: \code{"equi"} and \code{"unif"} for unperiodic sampling, \code{"sine"} and \code{"trian"} for sampling with a periodic density (see Details). 
}
  \item{npoints}{
integer: Sample size \eqn{n} (see Details).
}
  \item{ncycles}{
integer: Number of sampling cycles \eqn{n_s}{ns} (see Details).
}
  \item{ps}{
positive numeric value: Sampling period \eqn{p_s}{ps} (see Details).
}
}
\details{
\code{sampler} generates observation times \eqn{t_1,\ldots,t_n}{t} with a periodic sampling of period \eqn{p_s}{ps}. Four distributions are possible:
In case of \code{ttype="equi"}, the \eqn{t_i}{t[i]} are equidistantly sampled with \eqn{t_i=i\frac{p_sn_s}{n}}{t[i]=i*ps*ncycles/npoints}.
For \code{ttype="unif"}, the observation times are independently drawn form a uniform distribution on \eqn{[0,n_sp_s]}{[0,ncycles*ps]}.
Both these sampling schemes are aperiodic, the sampling period \eqn{p_s}{ps} only influences the length \eqn{t_n-t_1}{t[n]-t[1]} of the series of observation times.                       



For \code{ttype="sine"} and \code{ttype="trian"}, observation cycles \eqn{z^\star_i}{z.star[i]} are drawn from a uniform distribution on \eqn{\{1,\ldots,n_s\}}{{1,\ldots,ncycles}} and observation phases \eqn{\varphi^\star_i}{\phi.star[i]} are drawn from a density 
\deqn{d_{sine}(x)= \sin(2\pi x)+1}{d.sine(x)=sin(2\pi x)+1} (for \code{ttype="sine"}) or  
\deqn{d_{trian}(x)= 3x, \quad  0\leq x\leq\frac{2}{3},}{d.trian(x)=3x if 0 \le x \le 2/3} 
\deqn{d_{trian}(x)= 6-6x, \quad \frac{2}{3}<x\leq 1}{d.trian(x)=6-6x if 2/3 < x \le 1} (for \code{ttype="trian"}).
The unsorted observation times \eqn{t^\star_i}{t.star[i]} are then generated using 
\deqn{t^\star_i= \varphi^\star_i+(z^\star_i-1)p_s.}{t.star[i]=\phi.star[i]+(z.star[i]-1)*ps.}
Separately sampling observation cycle and phase was proposed by Hall and Yin (2003). For more details see Thieler, Fried and Rathjens (2016) or Thieler et al. (2013).
}
\value{
numeric vector: Ordered observation times.
}


\note{
To sample from \eqn{d_{sine}}{d.sine}, the function \code{BBsolve}, package \code{BB}, is used. 

A former version of this function is used in Thieler et al. (2013).
}


\author{
Anita M. Thieler and Jonathan Rathjens
}


\seealso{
Applied in \code{\link{tsgen}} (see there for an example).
}

\references{
Hall, P. and Yin, J. (2003): Nonparametric Methods for Deconvolving Multiperiodic Functions. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 65 (4), 869-886

Thieler, A. M., Backes, M., Fried, R. and Rhode, W. (2013): Periodicity Detection in Irregularly Sampled Light Curves by Robust Regression and Outlier Detection. Statistical Analysis and Data Mining, 6 (1), 73-89

Thieler, A. M., Fried, R. and Rathjens, J. (2016): RobPer: An R Package to Calculate Periodograms for Light Curves Based on Robust Regression. Journal of Statistical Software, 69 (9), 1-36, <doi:10.18637/jss.v069.i09>

}
