pcr.cv<-function(X,y,k=10,eps=0.000001,plot.it=FALSE){
    p<-ncol(X)
    cv <- rep(0,p+1)
    n <- nrow(X)
    all.folds <- split(sample(1:n), rep(1:k, length = n))
    coefficients.jackknife<-array(dim=c(p,p+1,k))
    for (i in seq(k)) {
        omit <- all.folds[[i]]
        Xtrain = X[-omit, , drop = FALSE]
        ytrain = y[-omit]
        Xtest = X[omit, , drop = FALSE]
        ytest = y[omit]
        pcr.object <- pcr(Xtrain,ytrain,eps)
        res <- matrix(, length(ytest),p+1)
        coefficients.jackknife[,,i]=pcr.object$coefficients
        for (j in 1:(p+1)){
            res[,j]<-ytest - pcr.object$intercept[j] - Xtest%*%pcr.object$coefficients[,j]
        
        }
        cv <- cv + apply(res^2, 2, sum)
    }
    cv <- cv/n
    m.opt <- which.min(cv)-1
    if (plot.it == TRUE) {
        plot(0:p, cv, type = "l")
    }
    pcr.object <- pcr(X,y,eps=eps)
    coefficients <- pcr.object$coefficients[,m.opt+1]
    intercept <- pcr.object$intercept[m.opt+1]
    cv.error<-min(cv)
    return(list(intercept = intercept, coefficients = coefficients, 
        m.opt = m.opt,cv.error=cv.error,coefficients.jackknife=coefficients.jackknife))




}
