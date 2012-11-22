pls.cv=function (X, y, k = 10, groups=NULL,m = ncol(X),use.kernel=FALSE,compute.covariance=FALSE) {
    n <- nrow(X)
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
    m = min(m, ntrain.min - 1,p)
    cv.error.matrix = matrix(0, k,m+1)
    rownames(cv.error.matrix)=my.names
    colnames(cv.error.matrix)=0:m
    for (i in seq(k)) {
        omit <- all.folds[[i]]
        Xtrain = X[-omit, , drop = FALSE]
        ytrain = y[-omit]
        Xtest = X[omit, , drop = FALSE]
        ytest = y[omit]
        pls.object <- pls.model(Xtrain, ytrain, m = m, Xtest = Xtest, 
            ytest = ytest, compute.DoF = FALSE, use.kernel=use.kernel)
        cv.error.matrix[i,] <-pls.object$mse
    }
    cv.error = apply(cv.error.matrix,2,mean)
    m.opt <- which.min(cv.error)-1
    if (compute.covariance==TRUE){
        use.kernel=FALSE
    }
    pls.object<-pls.model(X,y,m=max(m.opt,1),use.kernel=use.kernel,compute.DoF=compute.covariance,compute.jacobian=compute.covariance)
    intercept<-pls.object$intercept[m.opt+1]
    coefficients<-pls.object$coefficients[,m.opt+1]
    covariance<-pls.object$covariance
    if (compute.covariance==TRUE){
        covariance<-covariance[m.opt+1,,]
    }
    outlist=list(cv.error.matrix=cv.error.matrix,cv.error = cv.error, m.opt = m.opt,covariance=covariance,intercept=intercept,coefficients=coefficients)
    class(outlist) = "plsdof"
    return(outlist)
}
