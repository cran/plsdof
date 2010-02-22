gmdl=function (sigmahat, n, DoF, yhat) 
{
    #SS <- RSS/(n - DoF)
    SS<-sigmahat^2
    denominator<-DoF*SS
    FF <- (yhat)/(DoF * SS)
    FF[1,FF==0]=Inf
    
    dummy <- (n/2) * log(SS) + (DoF/2) * log(FF) + (1/2) * log(n)
    #par<- where.max(-dummy)
    par<-first.local.minimum(as.vector(dummy))
 
    return(list(score=dummy,par=par))
}
