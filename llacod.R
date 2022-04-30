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
chooseS <- function(s){
  0.01
}
cod_update <- function(w,theta0,y,X){
  n <- nrow(X)
  theta = theta0
  for (j in 1:ncol(X)) {
    r <- sum(t(y - X%*%theta+X[,j]*theta[j])%*%X[,j])/n
    z <- t(X[,j])%*%X[,j]/n
    g <- sign(r)*max(abs(r)-w[j],0)/z
    theta[j] <- g
  }
theta
}


cod <- function(w,beta0,y,X){
  theta0 <- beta0
  theta1 <- cod_update(w,theta0,y,X)
  count = 0
  print(n2(theta1-theta0))
  while((n2(theta1-theta0))>1e-5){
    theta0 <- theta1
    theta1 <- cod_update(w,theta0,y,X)
    count <- count+1
    print(n2(theta1-theta0))
    print(paste0("cod",count))
    if (count>5000) {
      print("fail cod")
      break
    }
  }
  if (count<5000) {
    print(paste0("success cod",count))
  }
  theta1
}





LLACOD <- function(beta.ini,lambda,a,y,X){
  beta0 = beta.ini
  w <- SCADdev(t = abs(beta0),lambda,a )
  beta1 <- cod(w,beta0,y,X)
  count = 0
  while (max(abs(beta1-beta0))>0.0001) {
    # print("newlla")
    beta0 <- beta1
    w <- SCADdev(t = abs(beta0),lambda,a )
    beta1 <- cod(w,beta0,y,X)
    count <- count+1
    
    print(paste0("***********LLA= ",count,"tol= ",max(abs(beta1-beta0))))
    if (count>5000) {
      print("fail LLA")
      break
    }
  }
  beta1
}


maincod <- function(data,
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
  LLACOD(beta.ini =beta0,lambda = lambda,a = a,y = y,X = X)
}

