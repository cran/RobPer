betaCvMfit <-
            # Robustly fits a Beta distribution to data using Cramer-von-Mises (CvM)
            # distance minimization.
function(data,CvM=TRUE, rob=TRUE){  
      if(any(is.na(data))) data <- data[-which(is.na(data))]  # delete NA-observations from data
      data[which(data<0)]<-0                                  # set negative observations to zero
      if(rob==FALSE){
            x.      <- mean(data)
            s..     <- var(data)}
      if(rob==TRUE){
            x.      <- median(data)
            s..     <- mad(data)^2
      }
      if(s..!=0){
      shape1  <- max(0.00001,-(x.*(-x.+x.^2+s..))/s..)  # moment estimators for the beta distribution forced to be positive
      shape2  <- max(0.00001,(shape1-shape1*x.)/x.)    
      erg<- c(shape1, shape2)
      if(CvM) {ff<-function(x) CvMbeta(x, vals=data)
                 shapes  <- optim(c(shape1, shape2), ff)$par  # minimize CvM distance
                 erg<- shapes}
      return(erg)
      }
      if(s..==0) return(c(NA, NA)) # if variance of the data is zero, no moment estimators can be calculated and no parameters are determined.
      }
