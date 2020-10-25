# Project Name
RICWL R package.

Robust index of confidence weighted learning for estimating optimal individual treatment rule when outcome has skwewed or heavy-tailed distributions or with outliers.
By using a similarity weighted indicator as our index of confidence, it uses only direction of paired outcomes but not their absoluate values, thus achived robustness towards irrgular outcome distributions.

## Table of contents
* [Setup](#setup)
* [Examples](#features)
* [Status](#status)
* [Contact](#contact)


## Setup
install.packages("devtools")
library(devtools)
install_github("sambiostat/RICWL")

## Code Examples
```r
library(RICWL)
library(caret)

#Generative Model 
sim<-function(N=500, p=5, sigma=0.3){
  X=matrix(runif(N*p,-1,1), nrow=N)
  colnames(X)=paste("X",c(1:p),sep="")
  A=sample(c(1,-1), N, replace=T, prob=c(0.5,0.5))
  Ym = 1+3*cos(X[,3])- 2*cos(X[,4])
  Yc = 4*(0.6- X[,1]^2 -X[,2]^2)

  #Homo cauchy
  noise=rcauchy(N, location = 0, scale = sigma)
 
  s2n=var(Yc) / (var(Ym)+ var(noise))
  Y=Ym + A*Yc + noise
  
  optA={Yc > 0}*2-1
  Y0=Ym - Yc #Everyone follows A=-1
  Y1=Ym + Yc #Everyone follows A=1
  optY=Ym + abs(Yc)
  
  out=list("X"=X, "A"=A, "Y"=Y, "optA"=optA,  "Y1"=Y1, "Y0"=Y0, "optY"=optY,
           "s2n"=s2n, "PrOptA1"= tmp[2]/N, "EY"=mean(optY), "EYopt0"=mean(optY[optA==-1]),
           "EYopt1"=mean(optY[optA==1]), "EYAll0"=mean(Y0), "EYAll1"=mean(Y1),
           "MedY"=median(optY), "MedYopt0"=median(optY[optA==-1]),
           "MedYopt1"=median(optY[optA==1]),"MedYAll0"=median(Y0),
           "MedYAll1"=median(Y1))

  return(out)
}

N=500
p=5
sgm=0.3
s=1
theta_list= c(0.25, 0.5, 0.75, 1)
sgm_list=c(0.001,0.01, 1,10,100,1000)
c_list=4^(-2:2)
m=5

test.data <- sim(N=10000, p=p, sigma=sgm, scenario = s)
X.test <- test.data$X
optTr.test <- factor(test.data$optA, levels=c(-1,1))
y.test <- test.data$y
Tr.test <- test.data$A
s2n.test <- test.data$s2n
y1.test <- test.data$Y1
y0.test <- test.data$Y0
optY.test <- test.data$EY
optYmed.test <- test.data$MedY

train.data <- sim(N=N, p=p, sigma=sgm, scenario = s)
X <- train.data[[1]]
Tr <- train.data[[2]]
y <- train.data[[3]]

#RICL_frac
pwl.frac.fit <- RICWL_CV(H=X, A=Tr, R2=y, pi=NULL, kernel='rbf', pentype = "lasso", XS=NULL, sigmalst=sgm_list,
                                   m=m,  e=1e-5, clist=c_list, theta_list = theta_list ,frac.par=0.5, method="frac")
  
optTr.pwl.frac<- predict(pwl.frac.fit, x=X.test)
optTr.pwl.frac <- factor(optTr.pwl.frac, levels=c(-1,1))
cfm.pwl.frac <- caret::confusionMatrix(as.factor(optTr.pwl.frac), reference=as.factor(optTr.test))
pcd.pwl.frac<- cfm.pwl.frac$overall["Accuracy"]  
value.pwl.frac <- mean(c(y1.test[optTr.pwl.frac==1], y0.test[optTr.pwl.frac==-1]))
med.pwl.frac <- median(c(y1.test[optTr.pwl.frac==1], y0.test[optTr.pwl.frac==-1]))   

  
#RICL_cos
pwl.cos.fit <- RICWL_CV(H=X, A=Tr, R2=y, pi=NULL, kernel='rbf', pentype = "lasso", XS=NULL, sigmalst=sgm_list,
                                        m=m,  e=1e-5, clist=c_list, theta_list = theta_list ,frac.par=0.5, method="cos")                            
  
optTr.pwl.cos<-predict(pwl.cos.fit, x=X.test)
optTr.pwl.cos <- factor(optTr.pwl.cos, levels=c(-1,1))
cfm.pwl.cos <- caret::confusionMatrix(as.factor(optTr.pwl.cos), reference=as.factor(optTr.test))
pcd.pwl.cos<- cfm.pwl.cos$overall["Accuracy"]
value.pwl.cos<- mean(c(y1.test[optTr.pwl.cos==1], y0.test[optTr.pwl.cos==-1]))
med.pwl.cos<- median(c(y1.test[optTr.pwl.cos==1], y0.test[optTr.pwl.cos==-1]))
  
##other rwl  
rwl.fit <- ROWL(X, Tr, y, pi=rep(0.5, N), kernel ="rbf", clinear=c_list, sigma=sgm_list,m=m, residual=T)
optTr.rwl<-  predict(rwl.fit, X.test)
optTr.rwl<- factor(optTr.rwl, levels=c(-1,1))
cfm.rwl <- caret::confusionMatrix(as.factor(optTr.rwl), reference=as.factor(optTr.test))
pcd.rwl<- cfm.rwl$overall["Accuracy"]
value.rwl <- mean(c(y1.test[optTr.rwl==1], y0.test[optTr.rwl==-1]))
med.rwl <- median(c(y1.test[optTr.rwl==1], y0.test[optTr.rwl==-1]))

#owl
owl.fit <- ROWL(X, Tr, y, pi=rep(0.5, N), kernel ="rbf", clinear=c_list, sigma=sgm_list,m=m, residual=F)
optTr.owl <- predict(owl.fit, X.test)
optTr.owl<-factor(optTr.owl, levels=c(-1,1))
cfm.owl <- caret::confusionMatrix(as.factor(optTr.owl), reference=as.factor(optTr.test))
pcd.owl<- cfm.owl$overall["Accuracy"]
value.owl <- mean(c(y1.test[optTr.owl==1], y0.test[optTr.owl==-1]))
med.owl <- median(c(y1.test[optTr.owl==1], y0.test[optTr.owl==-1]))   

  
   
##combine the res
result <-cbind( c(value.pwl.frac, pcd.pwl.frac, med.pwl.frac),
	     c(value.pwl.cos, pcd.pwl.cos, med.pwl.cos),
             c(value.owl, pcd.owl, med.owl),
             c(value.rwl, pcd.rwl, med.rwl),
             c(optY.test, 1, optYmed.test ))
  
 
result 
```

## Status
Project is: _finished_ 



## Contact
Created by jinchun zhang(jczhang0818@gmail.com) - feel free to contact me!