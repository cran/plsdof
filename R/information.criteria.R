information.criteria <-
function(RSS,DoF,yhat,sigmahat,n,DoF.max){
    steps=1
    # bound DoF by min(n-1,p)
    for (i in 1:nrow(DoF)){
            for (j in 1:ncol(DoF)){
                DoF[i,j]<-min(DoF.max,DoF[i,j])
                DoF[i,j]<-max(DoF[i,j],j)
            }
    }
    sigmahat.bak<-sigmahat
    sigmahat<-min(sigmahat)
    # information criteria based on our DoF estimate
    par.aic<-aic(RSS,n,DoF,sigmahat)
    par.bic<-bic(RSS,n,DoF,sigmahat)
    par.gmdl<-gmdl(RSS,n,DoF,yhat)
    #if (steps>1){
    #    for (j in 2:steps){
    #        par.aic<-aic(RSS,n,DoF,sigmahat.bak[par.aic[1],par.aic[2]])
    #        par.bic<-bic(RSS,n,DoF,sigmahat.bak[par.bic[1],par.bic[2]])
    #        }
    #}
    # naive estimates
    DoF.naive=matrix(1:ncol(RSS),nrow=ncol(RSS),ncol=nrow(RSS))
    DoF.naive=t(DoF.naive)
    sigma.naive.bak<-sqrt((RSS/(n-DoF.naive)))
    sigma.naive<-max(sigma.naive.bak)
    par.aic.naive<-aic(RSS,n,DoF.naive,sigma.naive)
    par.bic.naive<-bic(RSS,n,DoF.naive,sigma.naive)
    par.gmdl.naive<-gmdl(RSS,n,DoF.naive,yhat)
    if (steps>1){
        for (j in 1:steps){
        par.aic.naive<-aic(RSS,n,DoF.naive,sigma.naive.bak[par.aic.naive[1],par.aic.naive[2]])
        par.bic.naive<-bic(RSS,n,DoF.naive,sigma.naive.bak[par.bic.naive[1],par.bic.naive[2]])
        }
    }
    return(list(DoF=DoF,par.aic=par.aic,par.gmdl=par.gmdl,par.bic=par.bic,par.bic.naive=par.bic.naive,par.aic.naive=par.aic.naive,par.gmdl.naive=par.gmdl.naive))
}
