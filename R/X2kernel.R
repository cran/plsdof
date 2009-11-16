X2kernel <-
function(X,Xtest=NULL,sigma=1){
    s<-length(sigma)
    n<-nrow(X)
    Ktest=NULL
    K=array(dim=c(s,n,n))
    dummy.kernel<-X%*%t(X)
    v<-diag(dummy.kernel)
    dummy<- 2*dummy.kernel - v%*%t(rep(1,n)) - rep(1,n)%*%t(v)
    if (is.null(Xtest)==FALSE){
        ntest<-nrow(Xtest)
        Ktest=array(dim=c(s,ntest,n))
        vtest<-diag(Xtest%*%t(Xtest)) # this is probably slow
        dummy.test<-2* Xtest%*%t(X) - vtest%*%t(rep(1,n)) -rep(1,ntest)%*%t(v)
    }
    for (j in 1:s){
        K[j,,]<-exp((1/(sigma[j]^2))*dummy)
        if (is.null(Xtest)==FALSE){
            Ktest[j,,]=exp(((1/sigma[j]^2))*dummy.test)
        }
    }
    # center the data in feature space
    for (j in 1:s){
        Ktest.center=NULL
        v<-apply(K[j,,],2,mean)
        mean.v<-mean(v)
        inter<-rep(1,n)%*%t(v)
        K.center=K[j,,] - inter - t(inter) + matrix(mean.v,n,n)
        K[j,,]=K.center
        if (is.null(Ktest)==FALSE){
            row.mean<-apply(Ktest[j,,],1,mean)
            Ktest.center=Ktest[j,,]-row.mean%*%t(rep(1,n)) - rep(1,ntest)%*%t(v)+matrix(mean.v,ntest,n)
            Ktest[j,,]=Ktest.center
    }
    }
return(list(K=K,Ktest=Ktest,sigma=sigma))
}

