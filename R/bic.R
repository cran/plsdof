bic <-
function(RSS,n,DoF,sigmahat){
    dummy<-RSS/n + log(n)*(DoF/n)*sigmahat^2
    par.bic<-where.max(-dummy)
    return(par.bic)
}

