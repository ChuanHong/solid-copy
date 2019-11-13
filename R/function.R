g.logit=function (xx) 
{
  exp(xx)/(1 + exp(xx))
}

dg.logit = function(xx){exp(xx)/(exp(xx)+1)^2}

Est.ALASSO.Approx.GLMNET=function (ynew, xnew, bini, N.adj, lambda.grid, modBIC = T, N.inflate = N.adj) 
{
  w.b = 1/abs(bini)
  tmpfit = glmnet(x = xnew, y = ynew, family = "gaussian", 
                  penalty.factor = w.b, alpha = 1, lambda = lambda.grid, 
                  intercept = F)
  LL = apply((c(ynew) - predict(tmpfit, xnew, type = "response"))^2, 
             2, sum) * N.inflate
  if (modBIC) {
    BIC.lam = LL + min(N.adj^0.1, log(N.adj)) * tmpfit$df
    m.opt = which.min(BIC.lam)
    bhat.modBIC = tmpfit$beta[, m.opt]
    lamhat.modBIC = tmpfit$lambda[m.opt]
  }
  else {
    bhat.modBIC = lamhat.modBIC = NA
  }
  BIC.lam = LL + log(N.adj) * tmpfit$df
  m.opt = which.min(BIC.lam)
  bhat.BIC = tmpfit$beta[, m.opt]
  lamhat.BIC = tmpfit$lambda[m.opt]
  return(list(bhat.BIC = bhat.BIC, bhat.modBIC = bhat.modBIC, 
              lambda.BIC = lamhat.BIC, lambda.modBIC = lamhat.modBIC))
}

roc.fun=function(dat, betahat){
  yyi = g.logit(cbind(1,as.matrix(dat[,-1]))%*%betahat)
  roc = try(ROC.Est.FUN(dat[,1],yyi,yy0=0.5,seq(.01,.99,by=.01)),silent=T)
  roc
}

A.fun=function (bet, dat) 
{
  yy = dat[, 1]
  xx.vec = cbind(1, dat[, -1])
  -t(c(dg.logit(xx.vec %*% bet)) * xx.vec) %*% xx.vec/length(yy)
}

U.fun=function (bet, dat) 
{
  yy = dat[, 1]
  xx.vec = cbind(1, dat[, -1])
  c(t(c(yy - g.logit(xx.vec %*% bet))) %*% xx.vec)/length(yy)
}


iteration.ginv.fun=function (dat.list, bini, kk.list) {
  K = length(dat.list)
  Uini.list = sapply(kk.list, function(kk) {
    U.fun(bini, dat = dat.list[[kk]])
  })
  Aini.list = lapply(1:K, function(kk) {
    A.fun(bini, dat = dat.list[[kk]])
  })
  Ahat.ini = Reduce("+", Aini.list)/K
  bhat.list = -ginv(Ahat.ini) %*% Uini.list + bini
  list(b.k = bhat.list, Ahat = Ahat.ini)
}

ROC.Est.FUN <- function(Di,yyi,yy0,fpr0=NULL,wgti=NULL,yes.smooth=F)
{
  out.yy <- out.pp <- out.AUC <- out.TPR <- out.FPR <- out.PPV <- out.NPV <- NULL
  if(is.null(wgti)){wgti=rep(1,length(Di))}; yyi = as.matrix(yyi); pp=ncol(as.matrix(yyi));  
  mu0 = sum(wgti*(1-Di))/sum(wgti); mu1 = 1-mu0  
  for(k in 1:pp)
  {
    yy = yy0; 
    if(!is.null(fpr0)){
      tpr.all = S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth); 
      fpr.all = S.FUN(yyi[,k],Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth);
      TPR = approx(c(0,fpr.all,1),c(0,tpr.all,1),fpr0,method="linear",rule=2)$y; 
      TPR = c(S.FUN(yy0,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth), TPR); 
      yy = c(yy,Sinv.FUN(fpr0,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth))           
      FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)
    }else{
      TPR = S.FUN(yy,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth); 
      FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)
    }
    out.yy = cbind(out.yy, yy)
    out.pp = cbind(out.pp, S.FUN(yy,Yi=yyi[,k],wgti,yes.smooth=yes.smooth))
    out.TPR = cbind(out.TPR,  TPR);  out.FPR  <- cbind(out.FPR,  FPR)
    PPV <- 1/(1+FPR*mu0/(TPR*mu1)); NPV <- 1/(1+(1-TPR)*mu1/((1-FPR)*mu0))
    out.PPV <- cbind(out.PPV, PPV); out.NPV <- cbind(out.NPV, NPV)
    #AUC <- sum((sum.I(yyi[,k],"<=",Yi=yyi[,k],Vi=Di*wgti)+sum.I(yyi[,k],"<",Yi=yyi[,k],Vi=Di*wgti))*(1-Di)*wgti/2
    #             )/(sum((1-Di)*wgti)*sum(Di*wgti))
    AUC = sum(S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth)*(1-Di)*wgti)/sum((1-Di)*wgti)
    out.AUC <- c(out.AUC, AUC)
  }
  out = c(out.AUC,out.yy,out.pp,out.FPR,out.TPR,out.PPV,out.NPV)
  out
}

AUC.FUN = function(data)
{
  dd = data[,1]; xx = data[,2]; n0 = sum(1-dd); n1 = sum(dd) 
  x0 = xx[dd==0]; x1 = xx[dd==1]
  sum((sum.I(x0, "<=", x1)+sum.I(x0,"<",x1))/2)/(n0*n1)
}


S.FUN <- function(yy,Yi,Di,yes.smooth=F)
{
  if(yes.smooth){
    Y1i = Yi[Di==1]; n1 = sum(Di); bw = bw.nrd(Y1i)/n1^0.6
    c(t(rep(1/n1,n1))%*%pnorm((Y1i-VTM(yy,n1))/bw))
  }else{
    return((sum.I(yy,"<",Yi,Vi=Di)+sum.I(yy,"<=",Yi,Vi=Di))/sum(Di)/2)
  }
  ##sum.I(yy,"<=",Yi,Vi=Di)/sum(Di)
}

Sinv.FUN <- function(uu,Yi,Di,yes.smooth=F)
{
  yy0<-unique(sort(Yi,decreasing=T)); ss0 <- S.FUN(yy0,Yi,Di,yes.smooth=yes.smooth) 
  return(approx(ss0[!duplicated(ss0)],yy0[!duplicated(ss0)],uu,method="linear",f=0,rule=2)$y)
}

VTM<-function(vc, dm){
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}

sum.I <- function(yy,FUN,Yi,Vi=NULL)
  ## sum_i I(yy FUN Yi)Vi
  # Vi weight
{
  if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}
  # for each distinct ordered failure time t[j], number of Xi < t[j]
  pos <- rank(c(yy,Yi),ties.method='f')[1:length(yy)]-rank(yy,ties.method='f')    
  if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos # number of Xi>= t[j]
  if (!is.null(Vi)) {
    ## if FUN contains '=', tmpind is the order of decending
    if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else  tmpind <- order(Yi)
    ##Vi <- cumsum2(as.matrix(Vi)[tmpind,])
    Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
    return(rbind(0,Vi)[pos+1,])
  } else return(pos)
}

AUC = function(D,score,wgt=NULL){
  if(is.null(wgt)) {wgt=rep(1,length(D))}
  auc = sum(S.FUN(score,Yi=score,D*wgt,yes.smooth=F)*(1-D)*wgt)/sum((1-D)*wgt)
  return(auc)
}

ROC = function(D,score){
  roc = ROC.Est.FUN(D,score,0.5,seq(.01,.99,by=.01))
  roc = matrix(roc[-1],ncol=6)
  colnames(roc) = c("cutoff","est.pos.rate","FPR","TPR","PPV","NPV")
  return(roc)
}
