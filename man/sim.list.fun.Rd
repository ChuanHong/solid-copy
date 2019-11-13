\name{sim.list.fun}
\alias{sim.list.fun}
\title{Data generation}
\description{
Generate simulated data for SOLID method}
\usage{
sim.list.fun(nn,K,p.x,cor,beta0)
}
\arguments{
  \item{nn}{sample size}
  \item{K}{number of blocks}
  \item{p.x}{dimension of independent variables x}
  \item{cor}{correlation between x}
  \item{beta0}{beta coefficient}
}
\details{
This function generates the data for SOLID method} 

\value{
\code{dat.list} returns generated data
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
}
\keyword{data generation}
