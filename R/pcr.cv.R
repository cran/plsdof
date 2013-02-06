pcr.cv<-function(X,y,k=10,m=min(ncol(X),nrow(X)-1),groups=NULL,scale=TRUE,eps=0.000001,plot.it=FALSE,compute.jackknife=TRUE){
      n<-nrow(X)    
      p<-ncol(X)
   if (is.null(groups)==FALSE){
	f=as.factor(groups)
	k=length(levels(f))
	my.names=levels(f)
    }
    if (is.null(groups)==TRUE){
     f<-rep(1:k, length = n)
     my.names<-1:k
    }
    
    all.folds <- split(sample(1:n),f)
    ntrain = vector(length = k)
    for (i in 1:k) {
        ntrain[i] = n - length(all.folds[[i]])
    }
    ntrain.min = min(ntrain)
    m = min(m, ntrain.min - 1)
    cv.error.matrix = matrix(0, k,m+1)
    rownames(cv.error.matrix)=my.names
    colnames(cv.error.matrix)=0:m
    coefficients.jackknife=NULL
    if (compute.jackknife==TRUE){
    coefficients.jackknife<-array(dim=c(p,m+1,k))
    }
    for (i in seq(k)) {
        omit <- all.folds[[i]]
        Xtrain = X[-omit, , drop = FALSE]
        ytrain = y[-omit]
        Xtest = X[omit, , drop = FALSE]
        ytest = y[omit]
        pcr.object <- pcr(Xtrain,ytrain,scale,m,eps)
        res <- matrix(, length(ytest),m+1)
        if (compute.jackknife==TRUE){
        coefficients.jackknife[,,i]=pcr.object$coefficients
	}
        for (j in 1:(m+1)){
            res[,j]<-ytest - pcr.object$intercept[j] - Xtest%*%pcr.object$coefficients[,j]
        
        }
        cv.error.matrix[i,] <- apply(res^2, 2, sum)
    }
    cv.error <- apply(cv.error.matrix,2,mean)
    m.opt <- which.min(cv.error)-1
    if (plot.it == TRUE) {
        plot(0:m, cv.error, type = "l")
    }
    pcr.object <- pcr(X,y,scale,m.opt,eps=eps)
    coefficients <- pcr.object$coefficients[,m.opt+1]
    intercept <- pcr.object$intercept[m.opt+1]
    return(list(intercept = intercept, coefficients = coefficients, 
        m.opt = m.opt,cv.error.matrix=cv.error.matrix,cv.error=cv.error,coefficients.jackknife=coefficients.jackknife))




}
