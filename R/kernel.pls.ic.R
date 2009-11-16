kernel.pls.ic <-
function(X,y,m=ncol(X),type="vanilla",sigma=1,step.size=1){
        n<-nrow(X)
        DoF.max=n-1
        if (type=="vanilla"){
            DoF.max=min(DoF.max,ncol(X))
        }
        pls.object<-kernel.pls(X,y,m=m,type=type,compute.DoF=TRUE,sigma=sigma,step.size=step.size)
        RSS<-pls.object$RSS
        DoF<-pls.object$DoF
        yhat<-pls.object$yhat
        sigmahat<-(pls.object$sigmahat)
        # CHECK WHAT HAPPENS IF WE USE THE ABOVE DESCRIPTION
        #sigmahat=min(pls.object$sigmahat)
        if (type=="vanilla"){
            RSS=matrix(RSS,nrow=1)
            DoF=matrix(DoF,nrow=1)
            sigmahat=matrix(sigmahat,nrow=1)
            yhat=matrix(yhat,nrow=1)
        }
    ic<-information.criteria(RSS,DoF,yhat=yhat,sigmahat=sigmahat,n,DoF.max)
    DoF=ic$DoF
    ## compute optimal sigmas as well
    m.aic<-ic$par.aic[2]
    m.bic<-ic$par.bic[2]
    m.gmdl<-ic$par.gmdl[2]
    m.bic.naive<-ic$par.bic.naive[2]
    m.aic.naive<-ic$par.aic.naive[2]
    m.gmdl.naive<-ic$par.gmdl.naive[2]
    sigma.aic<-sigma[ic$par.aic[1]]
    sigma.bic<-sigma[ic$par.bic[1]]
    sigma.gmdl<-sigma[ic$par.gmdl[1]]
    sigma.bic.naive<-sigma[ic$par.bic.naive[1]]
    sigma.aic.naive<-sigma[ic$par.aic.naive[1]]
    sigma.gmdl.naive<-sigma[ic$par.gmdl.naive[1]]
    return(list(DoF=DoF,m.aic = m.aic,m.bic = m.bic,m.gmdl=m.gmdl,m.bic.naive = m.bic.naive,m.aic.naive=m.aic.naive,m.gmdl.naive=m.gmdl.naive,sigma.aic = sigma.aic,sigma.bic = sigma.bic,sigma.gmdl=sigma.gmdl,sigma.bic.naive = sigma.bic.naive,sigma.aic.naive=sigma.aic.naive,sigma.gmdl.naive=sigma.gmdl.naive))
}
