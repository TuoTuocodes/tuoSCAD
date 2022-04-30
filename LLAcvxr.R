getwd()
setwd("C:/work/r/ece720")
source("basic.R")
library(rqPen)
library(CVXR)
SCADdev <- function(t,lambda,a){
  coef2 <- max(a*lambda-t,0)/((a-1)*lambda)
  tem  <- ifelse(t<=lambda,1,coef2)
  lambda*tem
}

ISTA_cvxr <- function(lambda,beta0,y,X){
  y <- y
  X <- X
  w <- lambda
  l <- length(w)
  beta <- Variable(l)
  obj <- sum((y - X%*% beta)^2)/(2*length(y))+t(w)%*%abs(beta)
  prob <- Problem(Minimize(obj))
  result <- solve(prob)
  r = result$getValue(beta)
  return(as.numeric(r))
}



LLAcvxr <- function(beta.ini,lambda,a,y,X){
  beta0 = beta.ini
  w <- SCADdev(t = abs(beta0),lambda,a )
  beta1 <- ISTA_cvxr(w,beta0,y,X)
  count = 0
  while (max(abs(beta1-beta0))>1e-5) {
    # print("newlla")
    beta0 <- beta1
    w <- SCADdev(t = abs(beta0),lambda,a )
    beta1 <- ISTA_cvxr(w,beta0,y,X)
    count <- count+1
    print(paste0("***********LLA= ",count,"tol= ",max(abs(beta1-beta0))))
    if (count>5000) {
      print("fail LLA")
      break
    }
  }
  beta1
}






maincvxr <- function(data,
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
  LLAcvxr(beta.ini =beta0,lambda = lambda,a = a,y = y,X = X)
}


