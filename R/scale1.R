scale1 <-
function(u, b, cc, initial.sc=median(abs(u))/.6745) {
# subfunction for the Fast-S algorithm for linear regression originally published in 
# Salibian-Barrera, M. and Yohai, V. (2006): A Fast Algorithm for S-Regression Estimates.
# Journal of Computational and Graphical Statistics, 15 (2), 414-427
    # find the scale, full iterations
    max.it <- 200
    # magic number alert
    #sc <- median(abs(u))/.6745
    sc <- initial.sc
    i <- 0
    eps <- 1e-20
    # magic number alert
    err <- 1
    while( ( (i <- i+1) < max.it ) && (err > eps) ) {
        sc2 <- sqrt( sc^2 * mean( rho( u / sc, cc ) ) / b   )
        err <- abs(sc2/sc - 1)
        sc <- sc2
    }
    return(sc)
}
