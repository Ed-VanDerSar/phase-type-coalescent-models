##-------------------------------------------------------------
 library("partitions")
##-------------------------------------------------------------
##
## Purpose:
## This function finds the state space and calculates the rate matrix,
## from a number of linked sequences and a coalescent parameter.
##
## Name: RateMatAndLeftandStateSpace
## Authors: Hobolth, Siri-Jegousse, Bladt
## Date: October 2018
##
## Input:
## n: Number of linked sequences at time 0
## alpha: parameter of the Beta-coalescent
##
## Output:
## List consisting of
## RateMat: Rate matrix
## StateSpace: Matrix with rows corresponding to the states
## A state is a n dimensional vector (a1,...,an)
##
##----------------------------------------------------------------
RateMatandStateSpace <- function(n,alpha){
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
  ##----------------------------------------------------
  ## Intensity matrix
  ##----------------------------------------------------
  Rate<-matrix(0,ncol=dim,nrow=dim)
  ## Following thealgorithm 3.4.
  for (i in 2:dim){
    for (j in 1:(i-1)){
      ## establishing differences between two states
      c=Rmatrix[i,]-Rmatrix[j,]
      ## Identifyung if the two states are compatible
      sum1<-c%*%rep(1:n)
      check1<-sum1[1,1]
      ## check1==0 means that the size of the new blocks equals the size of the disappearing blocks
      ## Identifying how many new blocks are created
      w1<-ifelse(c>0,1,0)
      sum2<-c%*%w1
      check2<-sum2[1,1]
      ## check2==1 means that only one new block is created
      ## Disappearing blocks
      w2<-ifelse(c<0,1,0)
      cneg<--c*w2
      ##Fullfilling the rate matrix
      if (check1==0 & check2==1) {
        provrate<-beta(sum(cneg)-alpha,sum(Rmatrix[j,])-sum(cneg)+alpha)/beta(alpha,2-alpha)
        for (k in 1:n){
          provrate<-provrate*choose(Rmatrix[j,k],cneg[k])
        }
        Rate[j,i]<-provrate}
    }
  }
  ## Diagonal part of the matrix
  for (i in 1:dim){
    Rate[i,i]<--sum(Rate[i,])
  }
  return(list(RateMat=Rate, StateSpace=Rmatrix))
}
