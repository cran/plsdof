information.criteria=function (RSS, DoF, yhat, sigmahat, n,criterion="bic"){
    if (criterion=="aic"){
        dummy<-aic(RSS,n,DoF,sigmahat)
    }
    if (criterion=="bic"){
        dummy<-bic(RSS,n,DoF,sigmahat)
}
        
    if (criterion=="gmdl"){
        dummy<-gmdl(sigmahat, n, DoF, yhat)
}
    par <- dummy$par
    score<-dummy$score

    return(list(DoF = DoF, par = par, score=score))
}
