singlePerM <-
function(tt,yy,s,n,ii, pp,W,rho, regression, design, steps, var1, tol, genoudcontrol){
    # Used by RobPer 
    # Calculates one periodogram bar for given trial period pp
    # Input not checked as it is assumed that the function RobPer gives the right input


    
    # Initial Full model:
    # Build design matrix
    XX    <- Xgen(tt=tt,n=n,s=s, pp=pp, design=design, steps=steps)
    if(is.na(XX[1])){ return(NA)}else{
    
    # Initial LTS, small nsamp in order to be faster
    alpha<- ifelse(all(dim(XX)%%2==c(1,0)),0.5-1/(n-dim(XX)[2]-1),0.5)
    tempFULL <- try(ltsReg(yy~0+XX, mcd=FALSE, nsamp=50, alpha=alpha), silent=TRUE) #ltsReg seldomly crashes... 
          if(inherits(tempFULL, "try-error"))  tempFULL <- try(ltsReg(yy~0+XX, mcd=FALSE, nsamp=50, alpha=alpha), silent=TRUE) #give it a second chance
          if(inherits(tempFULL, "try-error"))  tempFULL <- try(ltsReg(yy~0+XX, mcd=FALSE, nsamp=50, alpha=alpha), silent=TRUE) # and a third time
          if(inherits(tempFULL, "try-error"))  tempFULL <- suppressWarnings(rq(yy~0+XX, tau=.5, method="br")) ## if still problems, use L1
                    
          rr <- as.vector(tempFULL$residuals)
      
    
    # Estimate of the standard deviation
    sigma <- ifelse(var1, 1, median(abs(rr[rr!=0]))/0.675)    
    
    
    # Initial L1 regression of the intercept model
    tempINT <- suppressWarnings(rq(yy~0+ii, tau=.5, method="br"))
         ee <- as.vector(tempINT$residuals)
    
    # IRWLS for intercept model
    gamma <- IRWLS(yy, matrix_=ii, W=W, residuals_=ee, scale_=sigma, tol=tol)
  
    
    # Full model continueing
    # Use genetic optimization in case of biweight regression in order to get a good starting point
    if(regression=="bisquare"){
        beta_ <- tempFULL$coeff 
        beta_gamma <- numeric(length(beta_))
        if(design %in% c("step", "stepB","splines")) beta_gamma<- beta_gamma+gamma
        if(design %in% c("sine", "fourier(2)", "fourier(3)")) beta_gamma[1]<- gamma   
        optf<- function(beta) sum(rho((yy-XX%*%beta)/sigma))
        beta_<-suppressWarnings(genoud(optf, nvars=length(beta_),starting.values=rbind(beta_, beta_gamma), pop.size=genoudcontrol$pop.size, max=FALSE, max.generations=genoudcontrol$max.generations, wait.generations=genoudcontrol$wait.generations, solution.tolerance=tol, print.level=0)$par)
        rr<- as.vector(yy-XX%*%beta_)
        }

    # IRWLS for full model       
    beta  <- IRWLS(yy, matrix_=XX, W=W, residuals_=rr, scale_=sigma, tol=tol)


    
    bar<- 1- sum(rho((yy-XX%*%beta)/sigma))/ sum(rho((yy-ii*gamma)/sigma))
    
    return(bar)
    }
    }
