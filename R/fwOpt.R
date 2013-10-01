fwOpt <-
function(x, cc)
{
# subfunction for the Fast-Tau algorithm for linear regression originally published in 
# Salibian-Barrera, M., Willems, G. und Zamar, R. (2008): The fast-tau estimator for regression.
# Journal of Computational and Graphical Statistics, 17 Nr. 3, 659-682.
  
	tmp <- (-1.944 / cc^2 + 1.728 * x^2 / cc^4 - 0.312 * x^4 / cc^6 + 0.016 * x^6 / cc^8) / 3.25
	tmp[abs(x) < 2*cc] <- 1 / (3.25*cc^2)
	tmp[abs(x) > 3*cc] <- 0
	tmp
	
}
