#simulation model
sim<-function(scenario=c(1,2,3,4,5), N=500, p=5, sigma=0.3){
  X=matrix(runif(N*p,-1,1), nrow=N)
  colnames(X)=paste("X",c(1:p),sep="")
  A=sample(c(1,-1), N, replace=T, prob=c(0.5,0.5))
  Ym = 1+3*cos(X[,3])- 2*cos(X[,4])
  Yc = 4*(0.6- X[,1]^2 -X[,2]^2)
  if(scenario==1){
    #Homo Normal
    noise=rnorm(N, mean=0, sd=sigma)
  }
  else if(scenario==2){
    #Homo cauchy
    noise=rcauchy(N, location = 0, scale = sigma)
  }
  else if(scenario==3){
    Ym= 1 + 3*cos(X[,3])-2*cos(X[,4]) + 40*(cos(X[,1])^50*cos(X[,2])^50) #SPWL is way better than the OWL
    Yc= 4*(0.6- X[,1]^2 -X[,2]^2)
    noise=rnorm(N, mean=0, sd=sigma)
  }
  else if(scenario==4){
    Ym= 1 + 3*cos(X[,3])-2*cos(X[,4]) + 20*(cos(X[,1])^100*cos(X[,2])^100)
    Yc= 10*(1-X[,1]^2 -X[,2]^2)*(X[,1]^2 + X[,2]^2 - 0.3)
    noise=rnorm(N, mean=0, sd=sigma)
  }
  else if(scenario==5){
    Ym= 1 + 3*cos(X[,3])-2*cos(X[,4]) + 60*(cos(X[,1])^100*cos(X[,2])^100)
    Yc= 10*(1-X[,1]^2 -X[,2]^2)*(X[,1]^2 + X[,2]^2 - 0.3)
    noise=rnorm(N, mean=0, sd=sigma)
  }else {break;}
  
  s2n=var(Yc) / (var(Ym)+ var(noise))
  Y=Ym + A*Yc + noise
  
  optA={Yc > 0}*2-1
  Y0=Ym - Yc #Everyone follows A=-1
  Y1=Ym + Yc #Everyone follows A=1
  optY=Ym + abs(Yc)
  
  tmp=table(optA)
  out=list("X"=X, "A"=A, "Y"=Y, "optA"=optA,  "Y1"=Y1, "Y0"=Y0, "optY"=optY,
           "s2n"=s2n, "PrOptA1"= tmp[2]/N, "EY"=mean(optY), "EYopt0"=mean(optY[optA==-1]),
           "EYopt1"=mean(optY[optA==1]), "EYAll0"=mean(Y0), "EYAll1"=mean(Y1),
           "MedY"=median(optY), "MedYopt0"=median(optY[optA==-1]),
           "MedYopt1"=median(optY[optA==1]),"MedYAll0"=median(Y0),
           "MedYAll1"=median(Y1))
  
  #return(unlist(out[8:17]))
  return(out)
}