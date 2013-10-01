f.w <-
# subfunction for the Fast-S algorithm for linear regression originally published in 
# Salibian-Barrera, M. and Yohai, V. (2006): A Fast Algorithm for S-Regression Estimates.
# Journal of Computational and Graphical Statistics, 15 (2), 414-427
function(u, cc)
{
    # weight function = psi(u)/u
    tmp <- (1 - (u/cc)^2)^2
    tmp[ abs(u/cc) > 1 ] <- 0
    return(tmp)
}
