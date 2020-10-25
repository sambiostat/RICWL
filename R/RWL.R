#' ROWL, a function for OWL and RWL.
#'
#' OWL and RWL are benchmark models in our simulation study
#' @param H: n by P feature matrix
#' @param A: treatment, takes value -1 and +1 with size n
#' @param R2: outcome or residual vector with length n
#' @param pi: propensity score with length n
#' @param pentype: penalty in the residual estimation process, useful or RWL only
#' @param kernel: kernel used in SVM
#' @param residual: True or False, if True then RWL is deployed, if FALSE then OWL is deployed. Default is True.
#' @param sigma: hyper-parameter in SVM rbf kernel
#' @param clinear: hyper-parameter in SVM rbf kernel
#' @param m: number of fold for cross validation in parameter tunning
#' @param e: least tolerated error
#' @import e1071
#' @import kernlab
#' @import glmnet
#' @export
#'
#' @examples
#' n = 200
#' p = 5
#' H = matrix(rnorm(n*p), nrow=n)
#' A = rbinom(n=n, size=1, prob=0.5)*2 - 1
#' R2 = 4*H[,1]^2 + A*(H[,2]^2 + H[,4]^2 < 1.5)
#' fit = ROWL(H, A, R2)
#' optimalITR_pred = predict(fit, H)
#'
ROWL<-function (H, A, R2, pi = rep(1, n), pentype = "lasso", kernel = "linear", residual=TRUE,
                sigma = c(0.03, 0.05, 0.07), clinear = 2^(-2:2), m = 4, e = 1e-05)
{
  npar = length(clinear)
  n = length(A)
  p = dim(H)[2]
  if (residual==TRUE & max(R2) != min(R2)) {
    if (pentype == "lasso") {
      cvfit = cv.glmnet(H, R2, nfolds = m)
      co = as.matrix(predict(cvfit, s = "lambda.min", type = "coeff"))
    }
    else if (pentype == "LSE") {
      co = coef(lm(R2 ~ H))
    }
    else stop(gettextf("'pentype' is the penalization type for the regression step of Olearning, the default is 'lasso',\nit can also be 'LSE' without penalization"))
    r = R2 - cbind(rep(1, n), H) %*% co
  }
  else r = R2
  rand = sample(m, n, replace = TRUE)
  r = r/pi

  if (kernel == "linear") {
    V = matrix(0, m, npar)
    for (i in 1:m) {
      this = (rand != i)
      X = H[this, ]
      Y = A[this]
      R = r[this]
      Xt = H[!this, ]
      Yt = A[!this]
      Rt = r[!this]
      for (j in 1:npar) {
        model = wsvm(X, Y, R, C = clinear[j], e = e)
        YP = predict(model, Xt)
        V[i, j] = sum(Rt * (YP == Yt))/sum(YP == Yt)
      }
    }
    mimi = colMeans(V)
    best = which.max(mimi)
    cbest = clinear[best]
    model = wsvm(H, A, r, C = cbest, e = e)
  }
  if (kernel == "rbf") {
    nsig = length(sigma)
    V = array(0, c(npar, nsig, m))
    for (i in 1:m) {
      this = (rand != i)
      X = H[this, ]
      Y = A[this]
      R = r[this]
      Xt = H[!this, ]
      Yt = A[!this]
      Rt = r[!this]
      for (j in 1:npar) {
        for (s in 1:nsig) {
          model = wsvm(X, Y, R, "rbf", sigma = sigma[s],
                       C = clinear[j], e = e)
          YP = predict(model, Xt)
          V[j, s, i] = sum(Rt * (YP == Yt))/sum(YP ==
                                                  Yt)
        }
      }
    }
    mimi = apply(V, c(1, 2), mean)
    best = which(mimi == max(mimi), arr.ind = TRUE)
    bestC = clinear[best[1]]
    bestSig = sigma[best[2]]
    print(bestC)
    print(bestSig)
    model = wsvm(H, A, r, "rbf", bestSig, C = bestC, e = e)
  }
  model
}

