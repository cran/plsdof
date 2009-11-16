where.max <-
function(M){
    dummy<-apply(M,1,max)
    row.max<-which.max(dummy)
    col.max<-which.max(M[row.max,])
    maxx<-c(row.max,col.max)
    return(maxx)

}

