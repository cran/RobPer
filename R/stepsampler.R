stepsampler <-
# subfunction for modified versions of the Fast-Tau algorithm and the Fast-S algorithm for linear regression
# Original References:
#
# Salibian-Barrera, M. and Yohai, V. (2006): A Fast Algorithm for S-Regression Estimates.
# Journal of Computational and Graphical Statistics, 15 (2), 414-427 
#
# Salibian-Barrera, M., Willems, G. und Zamar, R. (2008): The fast-tau estimator for regression.
# Journal of Computational and Graphical Statistics, 17 Nr. 3, 659-682.
#
# This subfunction is used in order to perform linear regression of step functions more efficiently.
function(hnr, nichtnull)
 sample(rep(nichtnull[which(nichtnull[,2]==hnr),1],2),1)
