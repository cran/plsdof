pcr<-function(X,y,eps=0.000001){
    p <- ncol(X)
    n <- nrow(X)
    m<-p
    Beta <- matrix(, p, m) # matrix of regression coefficients
    mean.y <- mean(y)
    y <- scale(y, scale = FALSE)
    mean.X <- apply(X, 2, mean)
    sd.X <- apply(X, 2, sd)
    sd.X[sd.X == 0] = 1
    X <- X - rep(1, nrow(X)) %*% t(mean.X)
    X <- X/(rep(1, nrow(X)) %*% t(sd.X))
    S<-t(X)%*%X
    b<-t(X)%*%y
    eig<-eigen(S)
    lambda<-eig$values
    U<-eig$vectors
    lambdainv<-1/lambda
    lambdainv[lambda<eps]=0
    #cat(paste("lambdainv is ",lambdainv,"\n"))
    dummy<-diag(lambdainv)%*%t(U)%*%b
    for (i in 1:m){
        Beta[,i]<-U[,1:i,drop=FALSE]%*%dummy[1:i]
    }
    coefficients <- matrix(0, p, m + 1)
    coefficients[, 2:(m + 1)] = Beta/(sd.X %*% t(rep(1, m)))
    intercept <- rep(mean.y, m + 1) - t(coefficients) %*% mean.X
    return(list(intercept=intercept,coefficients=coefficients))
}
