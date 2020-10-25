#' RICWL, robust confidence index weighted learning functions.
#'
#' This is confidence index weighted learning that takes patient's characteristic X,
#' treatment A and outcome Y and estimate an optimal individual treament rule.
#' @param H: n by P feature matrix
#' @param A: treatment, takes value -1 and +1 with size n
#' @param R2: outcome or residual vector with length n
#' @param pi: propensity score with length n, when RCT, it shoud be set as 0.5 for each subject
#' @param kernel: kernel used in SVM
#' @param pentype: penalty in the residual estimation process, useful or RWL only
#' @param XS: variable used in similarity calculation, if not provided, use all features
#' @param sigma: hyper-parameter in SVM rbf kernel
#' @param cpar: hyper-parameter in SVM rbf kernel
#' @param theta: hyper-parameter governing the definition of neighborhood. It's the quantile probability for similarity
#' @param frac.par: if "frac" similarity is used, then provide this fraction parameter, ranges from 0-1
#' @param method: either "cos" or "frac" for similarity definition
#' @param m: number of fold for cross validation in parameter tunning, default is 10
#' @param e: least tolerated error, default is 1e-5
#'
#' @examples
#' n = 200
#' p = 5
#' H = matrix(rnorm(n*p), nrow=n)
#' A = rbinom(n=n, size=1, prob=0.5)*2 - 1
#' R2 = 4*H[,1]^2 + A*(H[,2]^2 + H[,4]^2 < 1.5) + rcauchy(n)
#' fit = RICWL(H, A, R2)
#' optimalITR_pred = predict(fit, H)

RICWL <- function(H, A, R2, pi=rep(0.5,n), kernel='rbf', pentype = "lasso", XS=NULL,
                 sigma,  cpar, theta , frac.par=0.5, method="cos", m=10,  e=1e-5){
  if (is.null(pi)) pi= rep(1, length(A))

  #calculate the p_hat
  p.hat=Cal_Phat(x=H, y=R2, A=A, method=method, frac.par=frac.par, XS=XS, sim.par=theta)
  w=(p.hat - 0.5)/pi

  fitmodel = wsvm(as.matrix(H), A, w, "rbf", sigma =sigma , C = cpar)
  return(fitmodel)
}
