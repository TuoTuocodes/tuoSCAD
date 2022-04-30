library(lars)
library(glmnet)
getwd()
setwd("C:/work/r/ece720")
source("basic.R")
source("LLAcvxr.R")
source("llaista.R")
source("llaFista.R")
source("llacod.R")
library(rqPen)
library(CVXR)
data("diabetes")
y =as.vector(diabetes$y)
y= (y-mean(y))/var(y)
y1 =100*y
x1 = 100*diabetes$x
####linear regression
fit = lm(y1~x1)
round(coef(fit),2)
c  =fit$coefficients
####
library(BMA)
fitbma = iBMA.bicreg( x = as.data.frame(x1),
                    Y = y1, thresProbne0 = 5, verbose = TRUE, maxNvar = 30, nIter = 2)
summary(fitbma)
####LASSO
set.seed(2022)
cv = cv.glmnet(x =x1, y=y1, alpha = 1)
best.lambda = cv$lambda.min
best.lambda
fitLASSO = glmnet(x =x1, y=y1, alpha = 1,lambda = best.lambda)
round(coef(fitLASSO),3)
####LLA + CVXR
our = maincvxr(formula =y1~x1,lambda = fitSCAD$lambdaMin,a = 3.8)
round(our,3)
n2(our-fitSCAD$beta)
####LLA +ista
our2 = mainista(formula =y1~x1,lambda = fitSCAD$lambdaMin,a = 3.8)
round(our2,3)
n2(our2-fitSCAD$beta)
####LLA + FISTA
our3 = mainFista(formula =y1~x1,lambda = fitSCAD$lambdaMin,a = 3.8)
round(our3,3)
n2(our3-fitSCAD$beta)
####LLA + COD
our4 = maincod(formula =y1~x1,lambda = fitSCAD$lambdaMin,a = 3.8)
round(our4,3)
n2(our4-fitSCAD$beta)
#### ILAMM
library(ILAMM)
fitSCAD = cvNcvxReg(as.matrix(x1), y1, penalty = "SCAD")
as.vector(round( fitSCAD$beta,3))
####


