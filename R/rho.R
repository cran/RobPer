rho <-
# subfunction for the Fast-S algorithm for linear regression originally published in 
# Salibian-Barrera, M. and Yohai, V. (2006): A Fast Algorithm for S-Regression Estimates.
# Journal of Computational and Graphical Statistics, 15 (2), 414-427
function(u, cc=1.56) {
    w <- abs(u)<=cc
    v <- (u^2/(2)*(1-(u^2/(cc^2))+(u^4/(3*cc^4))))*w +(1-w)*(cc^2/6)
    v <- v*6/cc^2
    return(v)
}
