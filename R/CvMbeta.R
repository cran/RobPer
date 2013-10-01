CvMbeta <-
            # Minimization criterion to minimize the Cramer-von-Mises (CvM) distance between 
            # a Beta distribution with parameters (x[1], x[2]) and a data sample "vals".  
            # Based on R-Code provided by Brenton C. Clarke to minimize the CvM distance
            # between a Gamma distribution and a data sample.
            # See Clarke, McKinnon and Riley (2012): A fast robust method for fitting gamma
            # distributions. Statistical Papers, 53 Nr. 4, 1001-1014.
            # This function is applied by betaCvMfit.
function(x, vals){
            nn<- length(vals)
            return(sum((suppressWarnings(pbeta(vals, shape1=x[1], shape2=x[2]))-(rank(vals)-0.5)/nn)^2)/nn)
            #+ 1/(12*nn^2)
            }
