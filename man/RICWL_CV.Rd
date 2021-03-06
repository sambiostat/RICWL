% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RICWL_CV.R
\name{RICWL_CV}
\alias{RICWL_CV}
\title{RICWL_CV, robust confidence index weighted learning functions with corss-validation for hyper-parameter auto-tunning.}
\usage{
RICWL_CV(
  H,
  A,
  R2,
  pi = rep(0.5, n),
  kernel = "rbf",
  pentype = "lasso",
  XS = NULL,
  sigmalst = c(0.1, 1, 5),
  clist = 4^(-2:2),
  theta_list = c(0.25, 0.5, 0.75, 1),
  frac.par = 0.5,
  method = "cos",
  m = 10,
  e = 1e-05
)
}
\arguments{
\item{H:}{n by P feature matrix}

\item{A:}{treatment, takes value -1 and +1 with size n}

\item{R2:}{outcome or residual vector with length n}

\item{pi:}{propensity score with length n, when RCT, it shoud be set as 0.5 for each subject}

\item{kernel:}{kernel used in SVM}

\item{pentype:}{penalty in the residual estimation process, useful or RWL only}

\item{XS:}{variable used in similarity calculation, if not provided, use all features}

\item{sigmalst:}{candidate value for hyper-parameter sigma of rbf kernel}

\item{cpar:}{candidate value for hyper-parameter C of rbf kernel}

\item{theta_list:}{candidate value for hyper-parameter theta governing the definition of neighborhood}

\item{frac.par:}{if "frac" similarity is used, then provide this fraction parameter, ranges from 0-1}

\item{method:}{either "cos" or "frac" for similarity definition}

\item{m:}{number of fold for cross validation in parameter tunning, default is 10}

\item{e:}{least tolerated error, default is 1e-5}
}
\description{
This is the confidence index weighted learning functions that takes patient's characteristic X,
treatment A and outcome Y and estimate an optimal individual treament rule. Using cross-validation, the function will
choose the hyper-parameter in rbf kernel and neigborhood definition automatically.
}
\examples{
n = 200
p = 5
H = matrix(rnorm(n*p), nrow=n)
A = rbinom(n=n, size=1, prob=0.5)*2 - 1
R2 = 4*H[,1]^2 + A*(H[,2]^2 + H[,4]^2 < 1.5) + rcauchy(n)
fit = RICWL_CV(H, A, R2)
optimalITR_pred = predict(fit, H)

}
