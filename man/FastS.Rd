\encoding{latin1}
\name{FastS}
\alias{FastS}

\title{
S-Regression using the Fast-S-Algorithm
}
\description{
Performs S-Regression using the Fast-S-Algorithm.
}
\usage{
FastS(x, y, Scontrol=list(int = FALSE, N = 100, kk = 2, tt = 5, b= .5,
 cc = 1.547, seed=NULL), beta_gamma)
}

\arguments{
  \item{x}{
numeric \eqn{(n\times p)}{(n,p)}-matrix: Designmatrix.
}
  \item{y}{
numeric vector: \eqn{n} observations.
}
  \item{Scontrol}{
list of length seven: control parameters (see Details).
}
  \item{beta_gamma}{
numeric vector: Specifies one parameter candidate of length 
\eqn{p} (see Details). } } \details{ The Fast-S-Algorithm to 
efficiently perform S-Regression was published by 
Salibian-Barrera and Yohai (2006). It bases on starting with a set 
of \code{N} parameter candidates, locally optimizing them, but 
only with \code{kk} iterations, optimizing the \code{tt} 
best candidates to convergence and then choosing the best 
parameter candidate. The rho-function used is the biweight 
function with tuning parameter \code{cc}, the value \code{b} is 
set to the expected value of the rho-function applied to the 
residuals. The default \code{cc=1.547} and \code{b=.5} is 
chosen following Rousseeuw and Yohai (1984) to obtain an 
approximative breakdown point of 0.5. When setting \code{int} to \code{TRUE},
this adds an intercept column to the design matrix. For more details see 
Salibian-Barrera and Yohai (2006) or Thieler, Fried and Rathjens (2016). 

The R-function \code{FastS} used in \code{RobPer} is a slightly 
changed version of the R-code published in Salibian-Barrera and Yohai (2006). It was changed in order to work more efficiently, 
especially when fitting step functions, and to specify one 
parameter candidate in advance. For details see Thieler, Fried and Rathjens (2016). } 
\value{ \item{beta}{numeric vector: Fitted parameter vector.} 
\item{scale}{numeric: Value of the objective function} } 

\references{ 

Rousseeuw, P. J. and Yohai, V. J. (1984): Robust Regression by 
Means of S-estimators. In Franke, J., \enc{Härdle}{Haerdle}, W. und Martin, D. 
(eds.): Robust and Nonlinear Time Series Analysis. Berlin New 
York: Springer, Lecture Notes in Statistics No. 26, 256-272 

Salibian-Barrera, M. and Yohai, V. (2006): A Fast Algorithm for 
S-Regression Estimates. Journal of Computational and Graphical 
Statistics, 15 (2), 414-427 

Thieler, A. M., Fried, R. and Rathjens, J. (2016): RobPer: An R Package to Calculate Periodograms for Light Curves Based on Robust Regression. Journal of Statistical Software, 69 (9), 1-36, <doi:10.18637/jss.v069.i09>

}

\author{ Matias Salibian-Barrera and Victor Yohai, modified by 
Anita M. Thieler } 


\seealso{
Applied in \code{\link{RobPer}}. See \code{\link{FastTau}} for example.
}

