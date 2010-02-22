pls.cv=function (X, y, k = 10, m = ncol(X),use.kernel=FALSE,compute.covariance=FALSE) {
    n <- nrow(X)
    p<-ncol(X)
    all.folds <- split(sample(1:n), rep(1:k, length = n))
    ntrain = vector(length = k)
    for (i in 1:k) {
        ntrain[i] = n - length(all.folds[[i]])
    }
    ntrain.min = min(ntrain)
    m = min(m, ntrain.min - 1,p)
    cv.error = rep(0, m+1)
    for (i in seq(k)) {
        omit <- all.folds[[i]]
        Xtrain = X[-omit, , drop = FALSE]
        ytrain = y[-omit]
        Xtest = X[omit, , drop = FALSE]
        ytest = y[omit]
        pls.object <- pls.model(Xtrain, ytrain, m = m, Xtest = Xtest, 
            ytest = ytest, compute.DoF = FALSE, use.kernel=use.kernel)
        cv.error <- cv.error + (pls.object$mse) * length(ytest)
    }
    cv.error = cv.error/n
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
    outlist=list(cv.error = cv.error, m.opt = m.opt,covariance=covariance,intercept=intercept,coefficients=coefficients)
    class(outlist) = "plsdof"
    return(outlist)
}
