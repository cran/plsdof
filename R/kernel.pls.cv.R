kernel.pls.cv <-
function(X,y,k=10,m=ncol(X),type="vanilla",sigma=1){
        n<-nrow(X)
        all.folds <- split(sample(1:n), rep(1:k, length = n))
        ntrain = vector(length = k)
        for (i in 1:k) {
            ntrain[i] = n - length(all.folds[[i]])
        }
        ntrain.min = min(ntrain)
        m = min(m, ntrain.min - 1) # ensure that the number of components is not larger than ntrain
        cv.error=rep(0,m)
        if (type=="gaussian"){
            cv.error=matrix(0,length(sigma),m)
        }
        for (i in seq(k)) {
            omit <- all.folds[[i]]
            Xtrain = X[-omit, , drop = FALSE]
            ytrain = y[-omit]
            Xtest = X[omit, ,drop = FALSE]
            ytest = y[omit]
            pls.object<-kernel.pls(Xtrain,ytrain,m=m,Xtest=Xtest,ytest=ytest,compute.DoF=FALSE,step.size=1,sigma=sigma,type=type)
            cv.error<-cv.error+(pls.object$mse)*length(ytest)
        }
        cv.error=cv.error/n
        if (type=="vanilla"){
            m.opt<-which.min(cv.error)
            sigma.opt<-NULL
        }
        if (type=="gaussian"){
            dummy.cv<-apply(cv.error,2,min)
            m.opt<-which.min(dummy.cv)
            sigma.opt<-sigma[which.min(cv.error[,m.opt])]
        }
    return(list(cv.error=cv.error,m.opt=m.opt,sigma.opt=sigma.opt))
}

