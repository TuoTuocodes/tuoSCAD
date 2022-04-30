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
FISTA_update <- function(w,theta0,y0,t0,sk,y,X){
  s <- sk
  n <- nrow(X)
  res <- y - X%*%y0
  Z <- y0+(s/n)*t(X)%*%res
  signZ <- sign(Z)
  tem <- pmax(abs(Z)-s*w,0)
  #  print(tem)
  theta1 = signZ*tem
  t1 =0.5*(1+sqrt(1+4*t0^2))
  y1 = theta1 + (t0-1)/t1*(theta1-theta0)
  return(list(theta1 =theta1,y1=y1,t1=t1))
}


FISTA <- function(w,beta0,y,X){
  theta0 <- beta0
  t0 = 1
  y0 = theta0
  sk <- chooseS(s)
  FISTA = FISTA_update(w,theta0,y0,t0,sk,y,X)
  theta1 =FISTA$theta1
  y1 =FISTA$y1
  t1 =FISTA$t1
  count = 0
  print(n2(theta1-theta0))
  while((n2(theta1-theta0))>1e-5){
    theta0 <- theta1
    t0 = t1
    y0 = y1
    sk <- chooseS(s)
    FISTA = FISTA_update(w,theta0,y0,t0,sk,y,X)
    theta1 =FISTA$theta1
    y1 =FISTA$y1
    t1 =FISTA$t1
    count <- count+1
    print(paste0("FISTA",count))
    if (count>5000) {
      print("fail fista")
      break
    }
  }
  if (count<5000) {
    print(paste0("success fista",count))
  }
  theta1
}





LLAFista <- function(beta.ini,lambda,a,y,X){
  beta0 = beta.ini
  w <- SCADdev(t = abs(beta0),lambda,a )
  beta1 <- FISTA(w,beta0,y,X)
  count = 0
  while (max(abs(beta1-beta0))>1e-5) {
    # print("newlla")
    beta0 <- beta1
    w <- SCADdev(t = abs(beta0),lambda,a )
    beta1 <- FISTA(w,beta0,y,X)
    count <- count+1
    print(paste0("***********LLA= ",count,"tol= ",max(abs(beta1-beta0))))
    if (count>5000) {
      print("fail LLA")
      break
    }
  }
  beta1
}


mainFista <- function(data,
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
  LLAFista(beta.ini =beta0,lambda = lambda,a = a,y = y,X = X)
}

