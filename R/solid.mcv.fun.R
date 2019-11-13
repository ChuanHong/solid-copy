solid.mcv.fun=function (dat.list, niter, ridge=TRUE)
{ 
  lambda.grid=10^seq(-100,3,0.1)
  ##initial value
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
  ##screen
  update1=iteration.ginv.fun(dat.list, bini,1:K) 
  bnew1= apply(update1$b.k,1,mean)
  Ahalf1 = svd(-update1$Ahat); Ahalf1 = Ahalf1$u%*%diag(sqrt(Ahalf1$d))%*%t(Ahalf1$v)
  ynew1=Ahalf1%*%bnew1
  xnew1=Ahalf1
  bhat1.DCOS.ind=Est.ALASSO.Approx.GLMNET(ynew1,xnew1, bnew1,N.adj=N,lambda.grid)$bhat.modBIC
  id.var=which(bhat1.DCOS.ind!=0)
  dat.list.short = lapply(1:K,function(kk){dat.list[[kk]][,id.var]}) 
  bnew.update = bnew1[id.var]
  roc.mcv=betahat.mcv=NULL
  for(kk in 1:K){
    dat.list.t=dat.list.short[-kk]
    dat.v=dat.list.short[[kk]]
    bini.kk = coef(glmnet(dat.list.t[[1]][,-1],dat.list.t[[1]][,1],standardize=F,alpha=0,lambda=lambda.path,family="binomial"))[,length(lambda.path)]##need to give a sequence for lambda rather than single value##
    bnew.update=bini.kk
    for(ii in 1:(niter-1)){
      bnew=bnew.update
      update.DCOS.ind = iteration.ginv.fun(dat.list=dat.list.t,bini=bnew,kk.list=1:length(dat.list.t))
      bnew.update = apply(update.DCOS.ind$b.k,1,mean) 
    }
    Ahalf.DCOS.ind = svd(-update.DCOS.ind$Ahat); Ahalf.DCOS.ind = Ahalf.DCOS.ind$u%*%diag(sqrt(Ahalf.DCOS.ind$d))%*%t(Ahalf.DCOS.ind$v)
    betahat.t = Est.ALASSO.Approx.GLMNET(ynew=Ahalf.DCOS.ind%*%bnew.update,xnew=Ahalf.DCOS.ind,bini=bnew.update,N.adj=N,lambda.grid)$bhat.BIC
    betahat.mcv=rbind(betahat.mcv, betahat.t)
    roc.mcv=rbind(roc.mcv, roc.fun(dat.v, betahat.t))
  }
  betahat.mcv.long=matrix(0, nrow=K, ncol=dim(dat.list[[1]])[2])
  betahat.mcv.long[,id.var]=betahat.mcv
  colnames(betahat.mcv.long)=c("intercept", paste("b", c(1:(dim(betahat.mcv.long)[2]-1)), sep=""))
  auc.mcv=mean(roc.mcv[,1])
  auc.mcv.se=sd(roc.mcv[,1], na.rm=T)/sqrt(K)
  roc.mcv.mtx=matrix(colMeans(roc.mcv)[-1],ncol=6)[-1,]
  colnames(roc.mcv.mtx)=c("cut","p.pos","fpr","tpr","ppv","npv")
  return(list(auc.mcv=auc.mcv, auc.mcv.se=auc.mcv.se, roc.mcv.mtx=roc.mcv.mtx))
}  
