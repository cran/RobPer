loss.S <-
# subfunction for the Fast-S algorithm for linear regression originally published in 
# Salibian-Barrera, M. and Yohai, V. (2006): A Fast Algorithm for S-Regression Estimates.
# Journal of Computational and Graphical Statistics, 15 (2), 414-427
function(u,s,cc) mean(rho(u/s,cc) )
