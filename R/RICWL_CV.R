#' RICWL_CV, robust confidence index weighted learning functions with corss-validation for hyper-parameter auto-tunning.
#'
#' This is the confidence index weighted learning functions that takes patient's characteristic X,
#' treatment A and outcome Y and estimate an optimal individual treament rule. Using cross-validation, the function will
#' choose the hyper-parameter in rbf kernel and neigborhood definition automatically.
#'
#' @param H: n by P feature matrix
#' @param A: treatment, takes value -1 and +1 with size n
#' @param R2: outcome or residual vector with length n
#' @param pi: propensity score with length n, when RCT, it shoud be set as 0.5 for each subject
#' @param kernel: kernel used in SVM
#' @param pentype: penalty in the residual estimation process, useful or RWL only
#' @param XS: variable used in similarity calculation, if not provided, use all features
#' @param sigmalst: candidate value for hyper-parameter sigma of rbf kernel
#' @param cpar: candidate value for hyper-parameter C of rbf kernel
#' @param theta_list: candidate value for hyper-parameter theta governing the definition of neighborhood
#' @param frac.par: if "frac" similarity is used, then provide this fraction parameter, ranges from 0-1
#' @param method: either "cos" or "frac" for similarity definition
#' @param m: number of fold for cross validation in parameter tunning, default is 10
#' @param e: least tolerated error, default is 1e-5
#' @import kernlab
#' @import e1071
#'
#' @examples
#' n = 200
#' p = 5
#' H = matrix(rnorm(n*p), nrow=n)
#' A = rbinom(n=n, size=1, prob=0.5)*2 - 1
#' R2 = 4*H[,1]^2 + A*(H[,2]^2 + H[,4]^2 < 1.5) + rcauchy(n)
#' fit = RICWL_CV(H, A, R2)
#' optimalITR_pred = predict(fit, H)
#'

RICWL_CV <- function(H, A, R2, pi=rep(0.5,n), kernel='rbf', pentype = "lasso", XS=NULL,sigmalst=c(0.1, 1, 5),
                    clist=4.^(-2:2), theta_list = c(0.25, 0.5, 0.75, 1),frac.par=0.5, method="cos", m=10,  e=1e-5){
  n=length(A)
  p=dim(H)[2]
  if (is.null(pi)) pi=rep(0.5, n)

  # cross-validation
  rand = sample(m, n, replace = TRUE)
  if (kernel == "linear") {
    V = matrix(NA, m, npar)
    for (i in 1:m) {
      this = (rand != i)
      X = H[this, ]
      Y = A[this]
      W = w[this]
      Xt = H[!this, ]
      Yt = A[!this]
      Wt = w[!this]
      for (j in 1:npar) {
        test= try(wsvm(X, Y, W, C = clinear[j], e = e), silent=F)
        if (class(test) %in% "try-error" ){
          next
        }else{
          model= wsvm(X, Y, W, C = clinear[j], e = e)
          YP = predict(model, Xt)
          V[i, j] = sum(Wt * (YP == Yt))/sum(YP == Yt)
        }
      }
    }
    mimi = colMeans(V)
    best = which.max(mimi)
    cbest = ifelse(is.na(clinear[best]), 1, clinear[best])
    model = wsvm(H, A,w, C = cbest, e = e)
  }
  if (kernel == "rbf") {
    par_all = expand.grid(theta_list, sigmalst, clist)
    #print(dim(par_all))
    hypres=sapply( 1:nrow(par_all), function(jj){
      hpi = as.numeric(par_all[jj, ])
      tmpres=NULL
      for (i in 1:m){
        this = (rand != i)
        Xtrain = H[this, ]
        Rtrain = R2[this]
        Atrain = A[this]
        Pitrain = pi[this]

        Xtest = H[!this, ]
        Rtest = R2[!this]
        Atest = A[!this]
        pTest = Cal_Phat(x=Xtest, y=Rtest, A=Atest, method=method, frac.par=frac.par,  XS=XS, sim.par=hpi[1])
        wTest = (pTest - 0.5)/pi[!this]

        subFit= tryCatch(RICWL(H=Xtrain, A=Atrain, R2=Rtrain, pi=Pitrain, kernel=kernel, pentype = pentype, XS=XS,
                                   sigma=hpi[2], m=m,  e=1e-5, cpar=hpi[3], theta=hpi[1], frac.par=frac.par, method=method),
                         error=function(e) {e} )

        if ("error" %in% class(subFit)) {
          tmpres=c(tmpres,NA)
          print(subFit)
        }else{
          YP = predict(subFit, Xtest)
          #print("nonerror")
          tmpres = c(tmpres, sum(wTest * (YP == Atest))/sum(YP == Atest))
        } }
      return(mean(tmpres, na.rm=T))
    } )
    #print(paste("hypres", hypres))

    best = which(hypres == max(hypres, na.rm=T), arr.ind=T)
    bestTheta = ifelse( is.null(par_all[best, 1]), 0.5,par_all[best, 1])
    bestSigma = ifelse( is.null(par_all[best, 2]), 0.05, par_all[best, 2])
    bestC = ifelse( is.null(par_all[best, 3]), 1, par_all[best, 3])

    if (length(best)>1) best=best[1]
    print( par_all[best,] )

    model = RICWL(H=H, A=A, R2=R2, pi=pi, kernel=kernel, pentype = pentype, XS=XS,
                      sigma=bestSigma,  m=m,  e=1e-5, cpar=bestC, theta=bestTheta,
                      frac.par=frac.par, method=method)
  }
  model
}
