Mscale <-
function(u, b, c, initialsc) 
{
# subfunction for the Fast-Tau algorithm for linear regression originally published in 
# Salibian-Barrera, M., Willems, G. und Zamar, R. (2008): The fast-tau estimator for regression.
# Journal of Computational and Graphical Statistics, 17 Nr. 3, 659-682.
if (initialsc==0)
    initialsc = median(abs(u))/.6745
maxit <- 100
sc <- initialsc
i <- 0 
eps <- 1e-10
err <- 1
while  (( i < maxit ) & (err > eps)) {
    sc2 <- sqrt( sc^2 * mean(rhoOpt(u/sc,c)) / b)
    err <- abs(sc2/sc - 1)
    sc <- sc2
    i <- i+1
}

return(sc)

}
