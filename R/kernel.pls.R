kernel.pls <-
function(X,y,m=ncol(X),Xtest=NULL,ytest=NULL,compute.DoF=FALSE,type="vanilla",sigma=1,step.size=1){
    n<-nrow(X)
    DoF<-NULL
    mean.y<-mean(y)  
    y<-scale(y,scale=FALSE)
    mse<-NULL
    prediction<-NULL
    coefficients<-NULL
    sigmahat<-NULL
    RSS<-NULL
    intercept<-NULL
    ############################
    ###### vanilla kernel ######
    ############################
    if (type=="vanilla"){
        mean.X<-apply(X,2,mean)
        sd.X<-apply(X,2,sd)
        sd.X[sd.X==0]=1
        X<-X- rep(1,nrow(X))%*%t(mean.X)
        X<-X/(rep(1,nrow(X))%*%t(sd.X))
        #X<-scale(X)
        K<-X%*%t(X)   
        pls.object<-kernel.pls.fit(K,y,m,compute.DoF=compute.DoF,step.size=step.size)
        A<-pls.object$Alpha
        DoF<-pls.object$DoF
        coefficients<-t(X)%*%A
        coefficients = coefficients/(sd.X %*% t(rep(1, ncol(coefficients))))
        intercept<-rep(mean.y,m) -  t(coefficients) %*% mean.X
        Yhat<-pls.object$Yhat +mean.y
        yhat<-apply(Yhat^2,2,sum)
        # compute estimated sigma and RSS for the information criteria
        if (compute.DoF==TRUE){
            RSS<-pls.object$RSS
            sigmahat<-pls.object$sigmahat
        }
        # compute prediction for test data #
        if (is.null(Xtest)==FALSE){
            prediction= rep(1, nrow(Xtest)) %*% t(intercept) + Xtest%*%coefficients
            # compute mse on test data #
            if (is.null(ytest)==FALSE){
                res<-matrix(,nrow(Xtest),m)
                for (l in 1:m){
                    res[,l]=ytest-prediction[,l]
                }
                mse=apply(res^2,2,mean)
            }
        }
    }
    #############################
    ###### Gaussian kernel ######
    #############################
    if (type=="gaussian"){
        if (is.null(Xtest)==FALSE){
            prediction<-array(dim=c(length(sigma),nrow(Xtest),m))
            if (is.null(ytest)==FALSE){
                mse<-matrix(,length(sigma),m)
            }
        }
        K.test=NULL
        range.object<-myrange.X(X)
        X<-range.object$X
        if (is.null(Xtest)==FALSE){
            Xtest<-myrange.X(Xtest,range.object$a,range.object$b)$X
        }
        K.object<-X2kernel(X,Xtest,sigma=sigma)
        K<-K.object$K
        if (is.null(Xtest)==FALSE){
            Ktest<-K.object$Ktest
        }
        yhat<-matrix(,length(sigma),m)
        RSS<-yhat
        A<-array(dim=c(length(sigma),n,m))
        if (compute.DoF==TRUE){
            DoF<-matrix(,length(sigma),m)
            sigmahat<-matrix(,length(sigma),m)
            RSS<-matrix(,length(sigma),m)
        }
        for (s in 1:length(sigma)){
            pls.object<-kernel.pls.fit(K[s,,],y,m,compute.DoF=compute.DoF,step.size=step.size)
            A[s,,]=pls.object$Alpha
            Yhat<-pls.object$Yhat + mean.y
            yhat[s,]<-apply(Yhat^2,2,sum)
            if (compute.DoF==TRUE){
                sigmahat[s,]<-pls.object$sigmahat
                RSS[s,]<-pls.object$RSS
                DoF[s,]<-pls.object$DoF
            }
                # compute prediction for test data
               if (is.null(Xtest)==FALSE){
                prediction[s,,]= Ktest[s,,]%*%A[s,,] + mean.y
                # compute mse on test data #
                if (is.null(ytest)==FALSE){
                    res<-matrix(,nrow(Xtest),m)
                    for (l in 1:m){
                        res[,l]=ytest-prediction[s,,l]
                    }
                mse[s,]=apply(res^2,2,mean)
            }
        }
        }
    }
    
    
    return(list(prediction=prediction,mse=mse,coefficients=coefficients,intercept=intercept,A=A,DoF=DoF,RSS=RSS,Yhat=Yhat,sigmahat=sigmahat,yhat=yhat))
}

