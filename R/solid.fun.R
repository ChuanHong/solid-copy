solid.fun=function(dat.list, niter, ridge=TRUE)
{ 
  lambda.grid=10^seq(-100,3,0.1)
  K=length(dat.list)
  y1=dat.list[[1]][,1]
  x1=dat.list[[1]][,-1]
  p.x=dim(x1)[2]
  N=dim(do.call(rbind,dat.list))[1]
  lambda.ini=log(p.x)/sum(y1)
  ##initial
  if(ridge==F){
    bini = glm(y1~x1,family="binomial")$coef}
  if(ridge==T){
    lambda.path=seq(0.1,lambda.ini, length.out = 4)
    bini = coef(glmnet(x1,y1,standardize=F,alpha=0,lambda=lambda.path,family="binomial"))[,length(lambda.path)]##need to give a sequence for lambda rather than single value##
  }
  #update1
  update1=iteration.ginv.fun(dat.list=dat.list,bini=bini,kk.list=1:K)
  bnew1= apply(update1$b.k,1,mean)
  Ahalf1 = svd(-update1$Ahat); Ahalf1 = Ahalf1$u%*%diag(sqrt(Ahalf1$d))%*%t(Ahalf1$v)
  ynew1=Ahalf1%*%bnew1
  xnew1=Ahalf1
  bhat1.DCOS.ind=Est.ALASSO.Approx.GLMNET(ynew1,xnew1, bnew1,N.adj=N,lambda.grid)$bhat.modBIC 
  id.var=which(bhat1.DCOS.ind!=0)
  dat.list.short = lapply(1:K,function(kk){dat.list[[kk]][,id.var]}) 
  bnew.update = bnew1[id.var]
  
  for(ii in 1:(niter-1)){
    bnew=bnew.update
    update.DCOS.ind = iteration.ginv.fun(dat.list=dat.list.short,bini=bnew,kk.list=1:K)
    bnew.update = apply(update.DCOS.ind$b.k,1,mean) 
  }
  Ahalf.DCOS.ind = svd(-update.DCOS.ind$Ahat); Ahalf.DCOS.ind = Ahalf.DCOS.ind$u%*%diag(sqrt(Ahalf.DCOS.ind$d))%*%t(Ahalf.DCOS.ind$v)
  betahat.short = Est.ALASSO.Approx.GLMNET(ynew=Ahalf.DCOS.ind%*%bnew.update,xnew=Ahalf.DCOS.ind,bini=bnew.update,N.adj=N,lambda.grid)$bhat.BIC
  
  betahat.long=rep(0, dim(dat.list[[1]])[2])
  betahat.long[id.var]=betahat.short
  names(betahat.long)=c("intercept", paste("b", c(1:length(betahat.long[-1])), sep=""))
  res=list(betahat=betahat.long)
  res
}