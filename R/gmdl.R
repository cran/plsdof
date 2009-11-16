gmdl <-
function(RSS,n,DoF,yhat){
    SS<-RSS/(n-DoF)
    FF<- (yhat)/(DoF*SS)
    dummy<-(n/2)*log(SS) +(DoF/2)*log(FF) + (1/2)*log(n)
    par.gmdl<-where.max(-dummy)
    return(par.gmdl)
}

