\name{solid.fun}
\alias{solid.fun}
\title{solid function}
\description{
A screening and one-step linearization infused DAC}
\usage{
solid.fun(dat.list,niter, ridge)
}
\arguments{
  \item{dat.list}{dataset}
  \item{niter}{number of iterations}
  \item{ridge}{using ridge penalty for initial value}
}
\details{
Fit the SOLID method: a screening and one-step linearization infused DAC} 
\section{SOLID method}{Divide and conquer (DAC) is a commonly used strategy to overcome the challenges of extraordi- narily large data, by first breaking the dataset into series of data blocks, then combining results from individual data blocks to obtain a final estimation. Various DAC algorithms have been pro- posed to fit sparse predictive regression model in the L1 regularization setting. However, many existing DAC algorithms remain computationally intensive when sample size and number of can- didate predictors are both large. In addition, no existing DAC procedures provide inference for quantifying the accuracy of risk prediction models. In this paper, we propose a screening and one-step linearization infused DAC (SOLID) algorithm to fit sparse logistic regression to massive datasets, by integrating the DAC strategy with a screening step and sequences of linearization, which enables us to maximize the likelihood with only selected covariates and perform penalized estimation via a fast approximation to the likelihood. 
}
\value{
\code{betahat} returns estimated beta coefficients
}
\references{Hong, C., Wang, Y. and Cat T. (2019). A Divide-and-Conquer Method for Sparse Risk Prediction and Evaluation (under revision).
}
\author{Chuan Hong}

\examples{
N=1e5
p.x=50
K=10
n=N/K
niter=3
cor=0
b0=-8
bb = c(1, 0.8, 0.4, 0.2, 0.1)
beta0 = c(b0, bb, rep(0, p.x - length(bb)))
dat.list=sim.list.fun(nn=N,K=10,p.x=50,cor=0.2,beta0)
SOLID.fit=solid.fun(dat.list,niter=3)
}
\keyword{SOLID}
