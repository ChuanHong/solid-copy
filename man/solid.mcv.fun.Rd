\name{solid.mcv.fun}
\alias{solid.mcv.fun}
\title{Modified cross validation}
\description{
A screening and one-step linearization infused DAC}
\usage{
solid.mcv.fun(dat.list, niter, ridge)
}
\arguments{
  \item{dat.list}{dataset}
  \item{niter}{number of iterations}
  \item{ridge}{using ridge penalty for initial value}
}
\details{
Modified cross validation for SOLID method} 
\section{Modified Cross Validation}{To assess the accuracy of a predictive re- gression model, we develop a modified cross-validation (MCV) that utilizes the side products of the SOLID hence substantially reduce the computational burden. Compared with existing DAC methods, the MCV procedure is the first effort to make inference on accuracy.  
}
\value{
\code{auc.mcv} returns auc obtained by modified cross validation
\code{auc.mcv.se} returns the standard error for auc obtained by modified cross validation
\code{roc.mcv.mtx} returns the accuracy table obtained by modified cross validation
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
SOLID.mcv=solid.mcv.fun(dat.list, niter=3)
}
\keyword{modified cross validation}
