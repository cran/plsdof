dvvtz <-
function(v,z,dv,dz){
    if (is.matrix(v)==FALSE){
        v<-matrix(v,ncol=1)
        dv<-array(dv,dim=c(1,nrow(dv),ncol(dv)))
    }
    p=ncol(v)
    n<-nrow(v)
    dummy<-matrix(0,n,n)
    for (i in 1:p){
        D<-(v[,i]%*%t(z) + sum(v[,i]*z)*diag(n))%*%dv[i,,] + v[,i]%*%t(v[,i])%*%dz
        dummy<-dummy + D
    }
    return(dummy)
}
