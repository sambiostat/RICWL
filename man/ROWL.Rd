% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RWL.R
\name{ROWL}
\alias{ROWL}
\title{ROWL, a function for OWL and RWL.}
\usage{
ROWL(
  H,
  A,
  R2,
  pi = rep(1, n),
  pentype = "lasso",
  kernel = "linear",
  residual = TRUE,
  sigma = c(0.03, 0.05, 0.07),
  clinear = 2^(-2:2),
  m = 4,
  e = 1e-05
)
}
\arguments{
\item{H:}{n by P feature matrix}

\item{A:}{treatment, takes value -1 and +1 with size n}

\item{R2:}{outcome or residual vector with length n}

\item{pi:}{propensity score with length n}

\item{pentype:}{penalty in the residual estimation process, useful or RWL only}

\item{kernel:}{kernel used in SVM}

\item{residual:}{True or False, if True then RWL is deployed, if FALSE then OWL is deployed. Default is True.}

\item{sigma:}{hyper-parameter in SVM rbf kernel}

\item{clinear:}{hyper-parameter in SVM rbf kernel}

\item{m:}{number of fold for cross validation in parameter tunning}

\item{e:}{least tolerated error}
}
\description{
OWL and RWL are benchmark models in our simulation study
}
\examples{
n = 200
p = 5
H = matrix(rnorm(n*p), nrow=n)
A = rbinom(n=n, size=1, prob=0.5)*2 - 1
R2 = 4*H[,1]^2 + A*(H[,2]^2 + H[,4]^2 < 1.5)
fit = ROWL(H, A, R2)
optimalITR_pred = predict(fit, H)

}
