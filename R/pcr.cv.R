pcr.cv<-function(X,y,k=10,groups=NULL,scale=TRUE,eps=0.000001,plot.it=FALSE){
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
    cv.error.matrix = matrix(0, k,p+1)
    rownames(cv.error.matrix)=my.names
    colnames(cv.error.matrix)=0:p
    coefficients.jackknife<-array(dim=c(p,p+1,k))
    for (i in seq(k)) {
        omit <- all.folds[[i]]
        Xtrain = X[-omit, , drop = FALSE]
        ytrain = y[-omit]
        Xtest = X[omit, , drop = FALSE]
        ytest = y[omit]
        pcr.object <- pcr(Xtrain,ytrain,scale,eps)
        res <- matrix(, length(ytest),p+1)
        coefficients.jackknife[,,i]=pcr.object$coefficients
        for (j in 1:(p+1)){
            res[,j]<-ytest - pcr.object$intercept[j] - Xtest%*%pcr.object$coefficients[,j]
        
        }
        cv.error.matrix[i,] <- apply(res^2, 2, sum)
    }
    cv.error <- apply(cv.error.matrix,2,mean)
    m.opt <- which.min(cv.error)-1
    if (plot.it == TRUE) {
        plot(0:p, cv.error, type = "l")
    }
    pcr.object <- pcr(X,y,eps=eps)
    coefficients <- pcr.object$coefficients[,m.opt+1]
    intercept <- pcr.object$intercept[m.opt+1]
    return(list(intercept = intercept, coefficients = coefficients, 
        m.opt = m.opt,cv.error.matrix=cv.error.matrix,cv.error=cv.error,coefficients.jackknife=coefficients.jackknife))




}
