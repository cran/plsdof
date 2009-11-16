myrange.X <-
function(X,a=NULL,b=NULL){
    p<-ncol(X)
    if ((is.null(a)==TRUE) | (is.null(b)==TRUE)){
        a<-vector(length=p)
        b<-vector(length=p)
        min.X<-apply(X,2,min)
        max.X<-apply(X,2,max)
        for (i in 1:p){
            a[i]=-2/(min.X[i]-max.X[i])
            b[i]=1-a[i]*max.X[i]
        }
    }
    for (i in 1:p){
        X[,i]<-a[i]*X[,i]+b[i]
    }
    return(list(X=X,a=a,b=b))
}

