library(MASS)
##kernal
K <- function(u,type="Gaussian"){
  if(type == "Gaussian"){
    output <- pnorm(u,0,1)
  }
  return(output)
}
Kh <- function(u,type,h){
  K(u/h,type)
}
##sqrtM
sqrtM <- function(M){
  V <- eigen(M)$vectors
  D <- eigen(M)$values
  root <- V%*%diag(sqrt(D))%*%t(V)
}
#sqrtM(x)%*%sqrtM(x)
##  n2
n2 <- function(x){
  sqrt(sum(x*x))
}
#n2(c(1,0,1))
n1 <- function(x){
  max(abs(x))
}
## inter
dot <- function(x,y){
  sum(x*y)
}
##std
stdx <- function(x){
  mx <- x-mean(x)
  mx*sqrt(length(x))/(n2(x-mean(x)))
}
stdX <- function(X){
  apply(X,2,stdx)
}
##meanx
meanx <- function(x){
  x-mean(x)
}
meanX <- function(X){
  apply(X,2,meanx)
}