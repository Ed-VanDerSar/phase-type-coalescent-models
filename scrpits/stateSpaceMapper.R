library("partitions")

StateSpaceMapper <- function(n){
  ##----------------------------------------------------
  ## Possible states
  ##----------------------------------------------------
  ## Size of the state space
  dim<-P(n)
  ## Definition of the state matrix
  Rmatrix<-matrix(ncol=n,nrow=dim)
  ## Set of partitions of [n]
  x<-parts(n)
  ## Rewriting the partitions as (a1,...,an)
  for (i in 1:dim) {
    y<-x[,dim-i+1]
    for (j in 1:n){
      Rmatrix[i,j]<-length(which(y==j))
    }
  }
  ## Reordering
  Rmatrix<-Rmatrix[order(Rmatrix[,1],decreasing=TRUE),]
  
  return(Rmatrix)
}
