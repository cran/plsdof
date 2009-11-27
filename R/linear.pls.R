linear.pls<-function(X,y,m=ncol(X),model.selection="aic"){
    p<-ncol(X)
    n<-nrow(X)
    mean.X<-apply(X,2,mean)
    sd.X<-apply(X,2,sd)
    sd.X[sd.X==0]=1
    X<-X- rep(1,n)%*%t(mean.X)
    X<-X/(rep(1,n)%*%t(sd.X))
    mean.y<-mean(y)
    y<-scale(y,scale=FALSE)
    Beta<-matrix(,p,m)
    W<-V<-Beta
    dW<-dBeta<-dV<-array(dim=c(m,p,n))
    A<-t(X)%*%X
    b<-t(X)%*%y
    for (i in 1:m){
        if (i==1){
            W[,i]<-b
            dW[i,,]=t(X)
            dW[i,,]<-dA(W[,i],A,dW[i,,])
            dV[i,,]<-dW[i,,]
            W[,i]<-W[,i]/sqrt((sum((W[,i])*(A%*%W[,i]))))
            V[,i]<-W[,i]
            Beta[,i]<-sum(V[,i]*b)*V[,i]
            dBeta[i,,]<-dvvtz(V[,i],b,dV[i,,],t(X))
        }
        if (i>1){
            W[,i]<-b - A%*%Beta[,i-1]
            dW[i,,]<-t(X)- A%*%dBeta[i-1,,]
            V[,i]<-W[,i] - vvtz(V[,1:(i-1),drop=FALSE],A%*%W[,i])
            dV[i,,]=dW[i,,] - dvvtz(V[,1:(i-1),drop=FALSE],A%*%W[,i],dV[1:(i-1),,,drop=FALSE],A%*%dW[i,,])
            dV[i,,]<-dA(V[,i],A,dV[i,,])
            V[,i]<-V[,i]/sqrt((sum(t(V[,i])%*%A%*%V[,i])))
            Beta[,i]=Beta[,i-1] + sum(V[,i]*b)*V[,i]
            dBeta[i,,]<-dBeta[i-1,,]+ dvvtz(V[,i],b,dV[i,,],t(X))
        }
        }
        sigma<-RSS<-DoF<-yhat<-matrix(,1,m)
        covariance<-array(dim=c(m,p,p))
        DD<-diag(1/sd.X)
        for (i in 1:m){
            res<-y-X%*%Beta[,i]
            yhat[1,i]<-sum((X%*%Beta[,i])^2)
            RSS[1,i]<-sum(res^2)
            H<-X%*%dBeta[i,,]
            DoF[1,i]<-sum(diag(H))
            dummy<-(diag(n)-H)%*%(diag(n)-t(H))
            sigma[1,i]<-sqrt(RSS[1,i]/sum(diag(dummy)))
            covariance[i,,]<-sigma[1,i]^2*DD%*%dBeta[i,,]%*%t(dBeta[i,,])%*%DD
        }
        if (model.selection=="aic"){
            m.opt<-aic(RSS,n,DoF,sigma)[2]
        }
        if (model.selection=="bic"){
            m.opt<-bic(RSS,n,DoF,sigma)[2]
        }
        if (model.selection=="gmdl"){
            m.opt<-gmdl(RSS,n,DoF,yhat)[2]
        }
        coefficients<-Beta
        coefficients = coefficients/(sd.X %*% t(rep(1,m)))
        intercept<-rep(mean.y,m) -  t(coefficients) %*% mean.X 
        outlist=list(coefficients=coefficients,intercept=intercept,dBeta=dBeta,covariance=covariance,sigmahat=sigma,m.opt=m.opt,DoF=DoF)
        class(outlist)="plsdof"
    return(outlist)}