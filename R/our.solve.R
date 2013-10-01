our.solve <-
# subfunction for the Fast-S algorithm for linear regression originally published in 
# Salibian-Barrera, M. and Yohai, V. (2006): A Fast Algorithm for S-Regression Estimates.
# Journal of Computational and Graphical Statistics, 15 (2), 414-427
function(a,b) {
    a <- qr(a)
    da <- dim(a$qr)
    if(a$rank < (p <- da[2]))
        return(NA)
    else qr.coef(a, b)
}
