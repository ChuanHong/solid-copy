sim.list.fun=function(nn,K,p.x,cor,beta0)
{
  SIM.FUN=function (nn, p.x, cor,  beta0) 
  {
    Sig0.X = cor + (1 - cor) * diag(p.x)
    xx = mvrnorm(nn, mu = rep(0, p.x), Sigma = Sig0.X)
    prob.x = g.logit(c(cbind(1, xx) %*% beta0))
    yy = rbinom(nn, prob = prob.x, size = 1)
    cbind(yy, xx)
  }
  dat.list=list()
  n.k=nn/K
  for(k in 1:K){
    dat.list[[k]] = SIM.FUN(n.k, p.x = p.x, cor=cor, beta0)
  }
  dat.list
}
