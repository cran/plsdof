vcov.plsdof=function(object,...){
    return(object$covariance[object$m.opt,,])
}