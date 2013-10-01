WtellerOpt <-
function(x, cc)
{
# subfunction for the Fast-Tau algorithm for linear regression originally published in 
# Salibian-Barrera, M., Willems, G. und Zamar, R. (2008): The fast-tau estimator for regression.
# Journal of Computational and Graphical Statistics, 17 Nr. 3, 659-682.
	tmp <- (3.584 - 0.864 * x^4 / cc^4 + 0.208 * x^6 / cc^6 - 0.012 * x^8 / cc^8) / 3.25
	tmp[abs(x) < 2*cc] <- 0
	tmp[abs(x) > 3*cc] <- 2
	tmp
	
}
