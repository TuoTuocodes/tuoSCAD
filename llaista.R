getwd()
setwd("C:/work/r/ece720")
source("basic.R")
library(rqPen)
library(CVXR)
library(quantreg)
SCADdev <- function(t,lambda,a){
  coef2 <- max(a*lambda-t,0)/((a-1)*lambda)
  tem  <- ifelse(t<=lambda,1,coef2)
  lambda*tem
}
Q <- function(w,theta,y,X){
  obj <- sum((y - X%*% beta)^2)/(2*length(y))+t(w)%*%abs(beta)
  return(obj)
}
Qk <- function(w,theta,theta0,y,X){
  1
}
chooseS <- function(s){
  0.01
}
ISTA_update <- function(w,theta0,sk,y,X){
  s <- sk
  n <- nrow(X)
  res <- y - X%*%theta0
  Z <- theta0+(s/n)*t(X)%*%res
  signZ <- sign(Z)
  tem <- pmax(abs(Z)-s*w,0)
#  print(tem)
  signZ*tem
}


ISTA <- function(w,beta0,y,X){
  theta0 <- beta0
  sk <- chooseS(s)
  theta1 <- ISTA_update(w,theta0,sk,y,X)
  count = 0
  print(n2(theta1-theta0))
  while((n2(theta1-theta0))>1e-5){
    theta0 <- theta1
    sk <- chooseS(s)
    theta1 <- ISTA_update(w,theta0,sk,y,X)
    count <- count+1
    print(paste0("ISTA",count))
    if (count>5000) {
      print("fail ista")
      break
    }
  }
  if (count<5000) {
    print(paste0("success ista",count))
  }
  theta1
}





LLAista <- function(beta.ini,lambda,a,y,X){
  beta0 = beta.ini
  w <- SCADdev(t = abs(beta0),lambda,a )
  beta1 <- ISTA(w,beta0,y,X)
  count = 0
  while (max(abs(beta1-beta0))>1e-5) {
    # print("newlla")
    beta0 <- beta1
    w <- SCADdev(t = abs(beta0),lambda,a )
    beta1 <- ISTA(w,beta0,y,X)
    count <- count+1
    print(paste0("***********LLA= ",count,"tol= ",max(abs(beta1-beta0))))
    if (count>5000) {
      print("fail LLA")
      break
    }
  }
  beta1
}


mainista <- function(data,
                     formula,
                     a,
                     lambda,
                     T=5000)
{
  #####
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"),
             names(mf), 0)
  mf <- mf[c(1, m)]  ###### contain a space here 1
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval.parent(mf)
  mt <- attr(mf, "terms")
  X <- model.matrix(mt, mf)
  y <- model.extract(mf, "response")
  p <<- ncol(X)
  #  beta0 <- lm(y~X-1)$coefficients
  beta0 <- rep(0,p)
  ###
  LLAista(beta.ini =beta0,lambda = lambda,a = a,y = y,X = X)
}

