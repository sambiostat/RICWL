#' Cal_Phat:calculate the confidence index
#'
#' This is the first step in RICWL.
#' It calcualtes the confidence index based on the sum of similarity weighted indicator functions.
#' @param x:  n by P feature matrix
#' @param y: outcome or residual vector with sizse n
#' @param A: treatment, takes value -1 and +1 with size n
#' @param XS: names of feature should be used for similarity calculation, if NULL then use all
#' @param method: either cosine or fraction similarity, default is cosine "cos"
#' @param frac.par: the parameter needed for fraction similarity calculation
#' @param sim.par: hyper-parameter governing the definition of neighborhood. It is the minimal similairty quantile for neighborhood
#' @import lsa
#' @export
#' @examples
#' n = 200
#' p = 5
#' H = matrix(rnorm(n*p), nrow=n)
#' A = rbinom(n=n, size=1, prob=0.5)*2 - 1
#' R2 = 4*H[,1]^2 + A*(H[,2]^2 + H[,4]^2 < 1.5)
#' Cal_Phat(x=H, y=R2, A=A)

Cal_Phat<- function(x, y, A, XS=NULL, method=c("cos","frac"), frac.par=0.5, sim.par=0.2){

  n=length(A)

  if (!( method %in% c("cos","frac"))) { method="cos"}
  if(is.null(XS)) { H=as.matrix(x)} else{ H=as.matrix(x[, XS]) }

  if (method=="cos") {
    S=lsa::cosine( t(H) ) #input of cosine: row is the feature vectors, column is the instance;
    S.shift=(S+1)/2 #since cosine distance ranges from (-1) to (1).
    S=S.shift
  }

  if(method=="frac"){
    #fractional distance matrix is defined as
    #dist_d^f(x, y) = (\sum_i=1^d [ (x_i - y_i)^f] )^ (1/f)
    D=matrix(NA, nrow=n, ncol=n)
    for(i in 1:n){
      for (j in i:n){
        D[i,j]=D[j, i] = sum((abs(H[i, ] - H[j,]))^frac.par) ^ (1/frac.par)
      }}
    S=(max(D)-min(D))/(D)
  }

  #similarity is negatively proportional to the distance measurement,
  #print(paste("Method Used For the Similarity Measurement: ", method, sep=""))

  #Y could be the orignal outcome Y or the residual R.
  p.hat=rep(NA, n) #the probability of Y(a) > Y(-a)
  for (i in 1:n){
    Si.2qtl=quantile(S[i,], prob=sim.par) #cutoff for how much of the points will be included in the next step.
    B.idx=( A!=A[i] & S[i, ]> Si.2qtl)
    p.hat[i]=sum( ((S[i,])[B.idx]) * (y[i] > y[B.idx]))/sum((S[i,])[B.idx])
  }
  return(p.hat)
}
