pcr<-function (X, y, scale = TRUE,m=min(ncol(X),nrow(X)-1),eps = 1e-06)
{
    p <- ncol(X)
    n <- nrow(X)
    Beta <- matrix(, p, m)
    mean.y <- mean(y)
    #m=min(n,p)
    y <- scale(y, scale = FALSE)
    mean.X <- apply(X, 2, mean)
    if (scale == FALSE) {
        sd.X <- rep(1, p)
    }
    if (scale == TRUE) {
        sd.X <- apply(X, 2, sd)
        sd.X[sd.X == 0] = 1
    }
    X <- X - rep(1, nrow(X)) %*% t(mean.X)
    X <- X/(rep(1, nrow(X)) %*% t(sd.X))
    b <- t(X) %*% y
    	my.svd=svd(X)
    	lambda <- (my.svd$d)^2
    	U <- my.svd$v
    	lambdainv <- 1/lambda
    	lambdainv[lambda < eps] = 0
    	dummy <- diag(lambdainv) %*% t(U) %*% b
    	for (i in 1:m) {
		#cat(paste("component number ",i,"\n"))
        	Beta[, i] <- U[, 1:i, drop = FALSE] %*% dummy[1:i]
    	}
    	#if (p>m){
        #	for (i in ((m+1):(p))){
        #    		Beta[,i]=Beta[,m]
        #	}
    	#}
    
    coefficients <- matrix(0, p, m + 1)
    coefficients[, 2:(m + 1)] = Beta/(sd.X %*% t(rep(1, m)))
    intercept <- rep(mean.y, m + 1) - t(coefficients) %*% mean.X
    return(list(intercept = intercept, coefficients = coefficients))
}
