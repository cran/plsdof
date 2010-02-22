
bic=function(RSS, n, DoF, sigmahat) 
{
    dummy <- RSS/n + log(n) * (DoF/n) * sigmahat^2
    #par <- where.max(-dummy)
    par<-first.local.minimum(as.vector(dummy))
    return(list(score=dummy,par=par))
}
