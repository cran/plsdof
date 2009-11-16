aic <-
function(RSS,n,DoF,sigmahat){
    dummy<-RSS/n + 2*(DoF/n)*sigmahat^2
    par.aic<-where.max(-dummy)
    return(par.aic)
}

