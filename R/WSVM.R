
#' wsvm to solve the subject-weighted SVM
#' @param X: n by p feature matrix
#' @param A: label
#' @param wR: instance weight
#' @param kernel: propensity score with length n, can be either "rbf" or "linear", default is "linear"
#' @param sigma: hyper-parameter in SVM rbf kernel
#' @param C: hyper-parameter in SVM rbf kernel
#' @param e: least tolerated error
#' @import e1071
#' @import kernlab
#' @import glmnet
#' @return object in rbfcl or linearcf classs
#'


wsvm<-function(X, A, wR, kernel='linear',sigma=0.05,C=1,e=0.00001){
  wAR=A*wR
  if (kernel=='linear'){
    K=X%*%t(X)
  }
  else if (kernel=='rbf'){
    rbf=rbfdot(sigma=sigma)
    K=kernelMatrix(rbf,X)
  } else stop(gettextf("Kernel function should be 'linear' or 'rbf'"))

  H=K*(wAR%*%t(wAR))
  n=length(A)
  solution=ipop(-abs(wR),H,t(A*wR),0,numeric(n),rep(C,n),0,maxiter=100)
  alpha=primal(solution)
  alpha1=alpha*wR*A
  if (kernel=='linear'){
    w=t(X)%*%alpha1 #parameter for linear
    fitted=X%*%w
    rm=sign(wR)*A-fitted
  } else if (kernel=='rbf'){
    #there is no coefficient estimates for gaussian kernel
    #but there is fitted value, first we compute the fitted value without adjusting for bias
    fitted=K%*%alpha1
    rm=sign(wR)*A-fitted
  }
  Imid =(alpha < C-e) & (alpha > e)
  rmid=rm[Imid==1];
  if (sum(Imid)>0){
    bias=mean(rmid)
  } else {
    Iup=((alpha<e)&(A==-sign(wR)))|((alpha>C-e)&(A==sign(wR)))
    Ilow=((alpha<e)&(A==sign(wR)))|((alpha>C-e)&(A==-sign(wR)))
    rup=rm[Iup]
    rlow=rm[Ilow]
    bias=(min(rup)+max(rlow))/2}
  fit=bias+fitted
  if (kernel=='linear') {
    model=list(alpha1=alpha1,bias=bias,fit=fit,beta=w)
    class(model)<-'linearcl'
  } else if (kernel=='rbf') {
    model=list(alpha1=alpha1,bias=bias,fit=fit,sigma=sigma,X=X)
    class(model)<-'rbfcl'}
  return (model)
}


#' predict.linearcl: predict labels for a given new instance when model is trained using linear kernel
#' @param object: wsvm output under linear kernel
#' @param x: n by p feature matrix of new subjects
#' @return output label with length n
#' predict.linearcl()

predict.linearcl<-function(object,x,...){
  predict=sign(object$bias+x%*%object$beta)
  return(predict)
}


#' predict.rbfcl: predict labels for a given new instance when model is trained using linear kernel
#' @param object: wsvm output underl rbf kernel
#' @param x: n by p feature matrix of new subjects
#' @return output label with length n
predict.rbfcl<-function(object,x,...){
  rbf=rbfdot(sigma=object$sigma)
  n=dim(object$X)[1]
  if (is.matrix(x)) xm=dim(x)[1]
  else if (is.vector(x)) xm=1
  else stop('x must be vector or matrix')
  if (xm==1){ K <- apply(object$X,1,rbf,y=x)
  }else{   K<- matrix(0, xm, n)
  for (i in 1:xm) {K[i,]=apply(object$X,1,rbf,y=x[i,]) }}

  predict=sign(object$bias+K%*%object$alpha1)
  return(predict)
}
