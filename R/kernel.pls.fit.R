kernel.pls.fit <-
function(K,y,m,compute.DoF,step.size){
    # set variables
    n<-nrow(K)
    TT<-matrix(,n,m)
    Yhat<-matrix(,n,m)
    Alpha<-matrix(,n,m)
    Gamma<-matrix(,n,m)
    DoF=NULL
    dYhat=NULL
    if (compute.DoF==TRUE){
        dYhat=array(dim=c(m,n,n))
        dtildeT<-array(dim=c(m,n,n))
        dT<-dtildeT
        DoF=vector(length=m)
    }
    # compute all objects iteratively
    for (i in 1:m){
        if (i==1){
            ri<-y
            gi<-ri
            ti<-tildeti<-K%*%y
            if (compute.DoF==TRUE){
                dT[i,,]<-dtildeT[i,,]<-K
                dT[i,,]<-dnormalize(ti,dT[i,,])
                dtildeT[i,,]<-dT[i,,]
            }
            #ti<-normalize(ti)
            dummy<-normalize(ti,gi)
            ti<-dummy$v
            gi<-dummy$w
            tildeti<-ti
            TT[,i]<-ti
            Gamma[,i]<-gi
            Alpha[,i]<-Gamma[,i]*sum(ti*y)
            Yhat[,i]=vvtz(ti,y)
            if (compute.DoF==TRUE){
		ti<-as.vector(ti)
                dYhat[i,,]=dvvtz(ti,y,dT[i,,],diag(n))
                DoF[i]<-sum(diag(dYhat[i,,]))
            }
        }
        if (i>1){
            ri<-y-Yhat[,i-1]
            tildeti<-K%*%ri
            if (compute.DoF==TRUE){
                dtildeT[i,,]<-K%*%(diag(n)-dYhat[i-1,,])
            }
                j<-i-1
                if ((floor(i/step.size) == i/step.size)){
                   j=1 
                }
                TTi<-TT[,1:(i-1),drop=FALSE]
                ti<-tildeti- vvtz(TTi,tildeti)
                dummy<-rep(0,n)
                for (r in 1:(i-1)){
                    dummy<-dummy + sum(tildeti*TT[,r])*Gamma[,r]
                }
                gi<-ri- dummy
            if (compute.DoF==TRUE){
                    dT[i,,]=dtildeT[i,,]- dvvtz(TTi,tildeti,dT[1:(i-1),,,drop=FALSE],dtildeT[i,,])
                    dT[i,,]=dnormalize(ti,dT[i,,])
            }
            dummy<-normalize(ti,gi)
            ti<-dummy$v
            gi<-dummy$w
            #ti<-normalize(ti)
            TT[,i]<-ti
            Gamma[,i]<-gi
        Yhat[,i]<-Yhat[,i-1] + vvtz(ti,y)
        Alpha[,i]<-Alpha[,i-1] + Gamma[,i]*sum(TT[,i]*y)
        if (compute.DoF==TRUE){
	    ti<-as.vector(ti)
            dYhat[i,,]<- dYhat[i-1,,] + dvvtz(ti,y,dT[i,,],diag(n))
            DoF[i]<-sum(diag(dYhat[i,,]))
        }
    }
    
    }
    RSS=NULL
    sigmahat=NULL
    if ((compute.DoF)==TRUE){
        RES<-matrix(,n,m)
        for (i in 1:m){
        RES[,i]<-Yhat[,i]-y
        }
        RSS<-apply(RES^2,2,sum)
        denominator<-vector(length=m)
        for (i in 1:m){
            denominator[i]<-sum(diag((diag(n)-dYhat[i,,])%*%((diag(n)-t(dYhat[i,,])))))
            #DoF[i]=max(DoF[i],i)
            #DoF[i]=min(DoF[i],n-1)
            #denominator[i]=n-DoF[i]
        }
        sigmahat<-sqrt(RSS/denominator)
    }
    return(list(Alpha=Alpha,dYhat=dYhat,Yhat=Yhat,DoF=DoF,RSS=RSS,sigmahat=sigmahat))
}
