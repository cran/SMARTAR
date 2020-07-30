## ----setup--------------------------------------------------------------------
library(SMARTAR)
data(codiacs)

## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
seqmeans(data=codiacs ,family="gaussian",plot="d", digits = 2,xlab = "SMART design")
seqmeans(data=codiacs ,plot="s",color = "lightblue",xlab = "SEQ",family="gaussian")

## -----------------------------------------------------------------------------
atsmeans(data=codiacs,conf=TRUE, alpha=0.05,plot=TRUE,digits = 2,pch=18,xlab="Treatment sequence")
atsmeans(data=codiacs,conf=TRUE, alpha=0.05,digits = 2,pch=18,xlab="abc")

## -----------------------------------------------------------------------------
smartest(data=codiacs,method="IPW",adjust="Bon")

## -----------------------------------------------------------------------------
getncp(df=5, alpha = 0.05, beta = 0.2, d = 1e-04, start = 5)

## -----------------------------------------------------------------------------
smartsize(delta=0.0435,df=5,global=TRUE,alpha=0.05,beta=0.20)

SEQ <- 1:8
A1 <- c(rep(0,4),rep(1,4))
PI1 <- rep(0.5,8)
O2 <- rep(c(0,0,1,1),2)
P2 <- c(0.7,0.7,0.3,0.3,0.6,0.6,0.4,0.4)
A2 <- rep(c(0,1),4)
PI2 <- rep(0.5,8)
MEAN <- 1:8
SD <- rep(10,8)
SIMatrix <- as.data.frame(cbind(SEQ,A1,PI1,O2,P2,A2,PI2,MEAN,SD))
  
smartsize(SIMatrix,global=TRUE,alpha=0.05,beta=0.20)

