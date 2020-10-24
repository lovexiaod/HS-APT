###FDR test for correlation matrix in 4 different models
rm(list=ls(all=TRUE))
library(MASS);library(Matrix)
library(xtable)
library(plyr)
library(rgenoud)
library(fdrtool)
library(psych)
GI<-read.csv("G:/joint manuscript/GI.csv")
colnames(GI)[1] <- 'Month'
Mom<-read.csv("G:/joint manuscript/Mom.csv")
GI<-cbind(GI,rep(Mom[,2],each=11))
tdata<-dlply(GI,.(Month))
tdata<-tdata[(length(tdata)-137-24):length(tdata)]



nT=length(tdata)
N=11

Xlist<-list()
Ylist<-list()
for (i in 1:nT){
  Ylist[[i]]=tdata[[i]][,3]-tdata[[i]][,9]
  Xlist[[i]]=as.matrix(tdata[[i]][,c(4:8,10)])
}

#1MKT,2SMB,3HML,4RMW,5CMA,6MOM
c1<-c(1:5);c2<-c(1:3);c3<-c(1:3,6);c4<-c(1:4);c5<-c(1:4,6);c6<-c(1:6)
for (i in 1:nT){
  Xlist[[i]]=Xlist[[i]][,c1]
}
X1list<-list()
for (i in 1:nT){
  X1list[[i]]<-cbind(rep(1,N),Xlist[[i]])
}
p=ncol(Xlist[[1]])

Zlist<-list()
for (i in 1:nT){
  xx<-list()
  for (j in 1:N){
    xx[[j]]<-t(c(Xlist[[i]][j,]))
  }
  Zlist[[i]]<-as.matrix(bdiag(xx))
}

Z1list<-list()
for (i in 1:nT){
  xx<-list()
  for (j in 1:N){
    xx[[j]]<-t(c(1,Xlist[[i]][j,]))
  }
  Z1list[[i]]<-as.matrix(bdiag(xx))
}


W<-t(matrix(c(0,	0.08,	0.082,	0.116,	0.21,	0.083,	0.069,	0.127,	0.125,	0.049,	0.059,
              0.041,	0,	0.036,	0.236,	0.095,	0.026,	0.075,	0.05,	0.359,	0.036,	0.046,
              0.132,	0.113,	0,	0.147,	0.142,	0.067,	0.077,	0.08,	0.121,	0.057,	0.064,
              0.074,	0.295,	0.058,	0,	0.087,	0.031,	0.086,	0.064,	0.182,	0.052,	0.072,
              0.176,	0.157,	0.074,	0.114,	0,	0.051,	0.07,	0.079,	0.183,	0.043,	0.052,
              0.146,	0.09,	0.073,	0.086,	0.107,	0,	0.066,	0.198,	0.088,	0.067,	0.078,
              0.084,	0.179,	0.059,	0.164,	0.102,	0.046,	0,	0.074,	0.154,	0.063,	0.075,
              0.148,	0.114,	0.058,	0.117,	0.11,	0.131,	0.071,	0,	0.1,	0.066,	0.085,
              0.069,	0.388,	0.041,	0.157,	0.12,	0.027,	0.07,	0.047,	0,	0.035,	0.045,
              0.068,	0.097,	0.048,	0.113,	0.071,	0.052,	0.071,	0.078,	0.088,	0,	0.315,
              0.069,	0.107,	0.046,	0.132,	0.072,	0.052,	0.073,	0.085,	0.095,	0.268,	0),N,N))


##APT model with full covariance matrix 
XTX<-matrix(rep(0,N^2*p^2),N*p,N*p)
XTY<-matrix(rep(0,N*p),N*p,1)
XTX1<-matrix(rep(0,N^2*(p+1)^2),N*(p+1),N*(p+1))
XTY1<-matrix(rep(0,N*(p+1),N*(p+1),1))
for (i in 1:nT){
  XTX<-XTX+t(Zlist[[i]])%*%Zlist[[i]]
  XTY<-XTY+t(Zlist[[i]])%*%Ylist[[i]]
  XTX1<-XTX1+t(Z1list[[i]])%*%Z1list[[i]]
  XTY1<-XTY1+t(Z1list[[i]])%*%Ylist[[i]]
}
betahat<-solve(XTX)%*%XTY
betahat1<-solve(XTX1)%*%XTY1
sigma2hat<-0
sigma2hat1<-0
Sigmahat<-matrix(rep(0,N*N),N,N)
meanerr<-rep(0,N)
for (i in 1:nT){
  meanerr<-meanerr+Ylist[[i]]-Zlist[[i]]%*%betahat
}
meanerr<-meanerr/nT
yelist<-list()
for (i in 1:nT){
  sigma2hat<-sigma2hat+t(Ylist[[i]]-Zlist[[i]]%*%betahat)%*%(Ylist[[i]]-Zlist[[i]]%*%betahat)
  sigma2hat1<-sigma2hat1+t(Ylist[[i]]-Z1list[[i]]%*%betahat1)%*%(Ylist[[i]]-Z1list[[i]]%*%betahat1)
  Sigmahat<-Sigmahat+(Ylist[[i]]-Zlist[[i]]%*%betahat-meanerr)%*%t(Ylist[[i]]-Zlist[[i]]%*%betahat-meanerr)
  yelist[[i]]<-Ylist[[i]]-Zlist[[i]]%*%betahat
}
sigma2hat<-sigma2hat/(nT*N)
sigma2hat1<-sigma2hat1/(nT*N)
Sigmahat<-Sigmahat/nT
diagSigmahat<-diag(Sigmahat)
for (i in 1:N){
  for (j in 1:N){
    Sigmahat[i,j]<-Sigmahat[i,j]/sqrt(diagSigmahat[i]*diagSigmahat[j])
  }
}
Y_e<-data.frame(matrix(unlist(yelist), nrow=nT, byrow=T))

CE<-corr.test(Y_e)
CE$r[upper.tri(CE$r)]<-0
CE$p[upper.tri(CE$p)]<-1

pvalueSigma1<-CE$p[order(CE$p)][(N+1):((N+1)*N/2)]
fdr1<-p.adjust(pvalueSigma1)


fdrvector1<-rep(0,N*N)
fdrvector1[order(CE$p)[(N+1):((N+1)*N/2)]]<-fdr1
fdrmatrix1<-matrix(fdrvector1,N,N)

sum(round(fdrmatrix1,5)<0.01)-N-N*(N-1)/2
sum(round(fdrmatrix1,5)<0.05)-N-N*(N-1)/2

##S-APT model with homoscedasticity (kou 2018)
MWt1<-list();MBA11<-0;MBA12<-0;MBA21<-0;MBA22<-0;MBb1<-0;MBb2<-0;MBb3<-0
for (i in 1:nT){
  MWt1[[i]]<-W%*%Ylist[[i]]
  MBA11<-MBA11+t(MWt1[[i]])%*%MWt1[[i]]
  MBA12<-MBA12+t(MWt1[[i]])%*%Zlist[[i]]
  MBA21<-MBA21+t(Zlist[[i]])%*%MWt1[[i]]
  MBA22<-MBA22+t(Zlist[[i]])%*%Zlist[[i]]
  MBb1<-MBb1+t(MWt1[[i]])%*%Ylist[[i]]
  MBb2<-MBb2+t(Zlist[[i]])%*%Ylist[[i]]
  MBb3<-MBb3+t(Ylist[[i]])%*%Ylist[[i]]
}
LL_M<-function(rho){
  b<-solve(MBA22)%*%(MBb2-MBA21*rho)
  s<-(MBb3-2*rho*MBb1-2*t(b)%*%MBb2+rho*MBA11*rho+2*t(b)%*%MBA21*rho+t(b)%*%MBA22%*%b)/(N*nT)
  lc<--(N*nT)*log(2*pi*s)/2+nT*log(det((diag(1,N)-rho*t(W))%*%(diag(1,N)-W*rho)))/2-N*nT/2
  return(lc)
}

g_M<-genoud(LL_M, nvars=1,unif.seed=1234,max=TRUE, pop.size=1000,print.level = F)
rhohat2<-g_M$par
betahat2<-solve(MBA22)%*%(MBb2-MBA21*rhohat2)
yelist2<-list()
for (i in 1:nT){
  yelist2[[i]]<-Ylist[[i]]-rhohat2*W%*%Ylist[[i]]-Zlist[[i]]%*%betahat2
}
Y_e2<-data.frame(matrix(unlist(yelist2), nrow=nT, byrow=T))

CE2<-corr.test(Y_e2)
CE2$r[upper.tri(CE2$r)]<-0
CE2$p[upper.tri(CE2$p)]<-1


pvalueSigma2<-CE2$p[order(CE2$p)][(N+1):((N+1)*N/2)]
fdr2<-fdrtool(pvalueSigma2,statistic="pvalue")
fdrvector2<-rep(0,N*N)
fdrvector2[order(CE2$p)[(N+1):((N+1)*N/2)]]<-fdr2$qval
fdrmatrix2<-matrix(fdrvector2,N,N)
rownames(fdrmatrix2)<-c("Austria","Belgium","Finland","France","Germany","Greece","Ireland","Italy","Netherlands","Portugal","Spain")
colnames(fdrmatrix2)<-c("Austria","Belgium","Finland","France","Germany","Greece","Ireland","Italy","Netherlands","Portugal","Spain")
xtable(fdrmatrix2,digits=4)
sum(round(fdrmatrix2,5)<0.01)-N-N*(N-1)/2
sum(round(fdrmatrix2,5)<0.05)-N-N*(N-1)/2

##S-APT model with heteroscedasticity
Wt<-list()
for (l in 1:N){
  Wt[[l]]<-list()
}
A11<-list();A12<-list();A21<-list();A22<-list()
b1<-list()
b2<-list();b3<-list()
for (l in 1:N){
  A11[[l]]<-0;A12[[l]]<-0;A21[[l]]<-0;A22[[l]]<-0;b1[[l]]<-0;b2[[l]]<-0;b3[[l]]<-0
  for (i in 1:nT){
    Wt[[l]][[i]]<-W[l,]%*%Ylist[[i]] #1x1
    A11[[l]]<-A11[[l]]+(Wt[[l]][[i]])^2 #1x1
    A21[[l]]<-A21[[l]]+Xlist[[i]][l,]%*%Wt[[l]][[i]] #Px1
    A22[[l]]<-A22[[l]]+Xlist[[i]][l,]%*%t(Xlist[[i]][l,]) #PxP same
    b1[[l]]<-b1[[l]]+Ylist[[i]][l]*Wt[[l]][[i]] #1x1
    b2[[l]]<-b2[[l]]+Xlist[[i]][l,]*Ylist[[i]][l] #Px1 same
    b3[[l]]<-b3[[l]]+Ylist[[i]][[l]]*Ylist[[i]][[l]] #1x1 same
  }
}


DH_H2<-function(rho){
  b<-list()
  for (l in 1:N){
    b[[l]]<-solve(A22[[l]])%*%(b2[[l]]-A21[[l]]*rho) #Px1
  }
  s<-rep(0,N)
  for (l in 1:N){
    s[l]<-(b3[[l]]-2*b1[[l]]*rho-2*t(b[[l]])%*%b2[[l]]+(rho^2)*A11[[l]]+2*t(b[[l]])%*%A21[[l]]*rho+t(b[[l]])%*%A22[[l]]%*%b[[l]])/(nT)
  }
  lc<-(-N*nT*log(2*pi))/2-nT*log(det(diag(s)))/2+nT*log(det((diag(1,N)-t(W)%*%t(diag(rep(rho,N))))%*%(diag(1,N)-diag(rep(rho,N))%*%W)))/2-N*nT/2
  #lc<-(-N*nT*log(2*pi))/2-nT*sum(log(det(diag(s))))/2+nT*log(det((diag(1,N)-t(diag(rep(rho,N)))%*%t(W))%*%(diag(1,N)-W%*%diag(rep(rho,N)))))/2-N*nT/2
  return(lc)
}


parlist2<-list()
betalist2<-list()
likelihoodlist2<-list()
threepartsofL2<-list()
nrep2<-10
#set.seed(1234)
for (i in 1:nrep2){
  print(c(i,i))
  g_H<-genoud(DH_H2, nvars=1, max=TRUE, pop.size=1000,print.level = F)
  likelihoodlist2[[i]]<-g_H$value
  parlist2[[i]]<-g_H$par
  #threepartsofL[[i]]<-c(-(N*nT)*log(slist[[i]])/2,nT*log(det((diag(1,N)-t(diag(c(parlist[[i]])))%*%t(W))%*%(diag(1,N)-W%*%diag(c(parlist[[i]])))))/2,-N*nT/2)
  likelihoodlist2[[i]]<-g_H$value
}
rhohat3<-parlist2[[which.max(likelihoodlist2)]]
likeli<-likelihoodlist2[[which.max(likelihoodlist2)]]

b22<-list()
for (l in 1:N){
  b22[[l]]<-solve(A22[[l]])%*%(b2[[l]]-A21[[l]]*rhohat3) 
}

betahat3<-rep(0,N*p)
for (l in 1:N){
  betahat3[(1+(l-1)*p):(l*p)]<-b22[[l]]
}

yelist3<-list()
for (i in 1:nT){
  yelist3[[i]]<-Ylist[[i]]-rhohat3*W%*%Ylist[[i]]-Zlist[[i]]%*%betahat3
}
Y_e3<-data.frame(matrix(unlist(yelist3), nrow=nT, byrow=T))

CE3<-corr.test(Y_e3)
CE3$r[upper.tri(CE3$r)]<-0
CE3$p[upper.tri(CE3$p)]<-1

pvalueSigma3<-CE3$p[order(CE3$p)][(N+1):((N+1)*N/2)]
fdr3<-fdrtool(pvalueSigma3,statistic="pvalue")
sum(round(fdr3$qval,5)<0.01)
sum(round(fdr3$qval,5)<0.05)


##HS-APT model with heteroscedasticity
Wt<-list()
for (l in 1:N){
  Wt[[l]]<-list()
}
A11<-list();A12<-list();A21<-list();A22<-list()
b1<-list()
b2<-list();b3<-list()
for (l in 1:N){
  A11[[l]]<-0;A12[[l]]<-0;A21[[l]]<-0;A22[[l]]<-0;b1[[l]]<-0;b2[[l]]<-0;b3[[l]]<-0
  for (i in 1:nT){
    Wt[[l]][[i]]<-W[l,]%*%Ylist[[i]] #1x1
    A11[[l]]<-A11[[l]]+(Wt[[l]][[i]])^2 #1x1
    A21[[l]]<-A21[[l]]+Xlist[[i]][l,]%*%Wt[[l]][[i]] #Px1
    A22[[l]]<-A22[[l]]+Xlist[[i]][l,]%*%t(Xlist[[i]][l,]) #PxP same
    b1[[l]]<-b1[[l]]+Ylist[[i]][l]*Wt[[l]][[i]] #1x1
    b2[[l]]<-b2[[l]]+Xlist[[i]][l,]*Ylist[[i]][l] #Px1 same
    b3[[l]]<-b3[[l]]+Ylist[[i]][[l]]*Ylist[[i]][[l]] #1x1 same
  }
}


DH_H2<-function(rho){
  b<-list()
  for (l in 1:N){
    b[[l]]<-solve(A22[[l]])%*%(b2[[l]]-A21[[l]]*rho[l]) #Px1
  }
  s<-rep(0,N)
  for (l in 1:N){
    s[l]<-(b3[[l]]-2*b1[[l]]*rho[l]-2*t(b[[l]])%*%b2[[l]]+(rho[l]^2)*A11[[l]]+2*t(b[[l]])%*%A21[[l]]*rho[l]+t(b[[l]])%*%A22[[l]]%*%b[[l]])/(nT)
  }
  lc<-(-N*nT*log(2*pi))/2-nT*log(det(diag(s)))/2+nT*log(det((diag(1,N)-t(W)%*%t(diag(c(rho))))%*%(diag(1,N)-diag(c(rho))%*%W)))/2-N*nT/2
  #lc<-(-N*nT*log(2*pi))/2-nT*sum(log(det(diag(s))))/2+nT*log(det((diag(1,N)-t(diag(rep(rho,N)))%*%t(W))%*%(diag(1,N)-W%*%diag(rep(rho,N)))))/2-N*nT/2
  return(lc)
}


parlist2<-list()
betalist2<-list()
likelihoodlist2<-list()
threepartsofL2<-list()
nrep2<-10
#set.seed(1234)
for (i in 1:nrep2){
  print(c(i,i))
  g_H<-genoud(DH_H2, nvars=N, max=TRUE, pop.size=1000,print.level = F)
  likelihoodlist2[[i]]<-g_H$value
  parlist2[[i]]<-g_H$par
  #threepartsofL[[i]]<-c(-(N*nT)*log(slist[[i]])/2,nT*log(det((diag(1,N)-t(diag(c(parlist[[i]])))%*%t(W))%*%(diag(1,N)-W%*%diag(c(parlist[[i]])))))/2,-N*nT/2)
  likelihoodlist2[[i]]<-g_H$value
}
rhohat<-parlist2[[which.max(likelihoodlist2)]]

b22<-list()
for (l in 1:N){
  b22[[l]]<-solve(A22[[l]])%*%(b2[[l]]-A21[[l]]*rhohat[l]) 
}

betahat<-rep(0,N*p)
for (l in 1:N){
  betahat[(1+(l-1)*p):(l*p)]<-b22[[l]]
}

s<-rep(0,N)
for (l in 1:N){
  s[l]<-(b3[[l]]-2*t(rhohat[l])%*%b1[[l]]-2*t(b22[[l]])%*%b2[[l]]+t(rhohat[l])%*%A11[[l]]%*%rhohat[l]+2*t(b22[[l]])%*%A21[[l]]%*%rhohat[l]+t(b22[[l]])%*%A22[[l]]%*%b22[[l]])/(nT)
}

Y_e2<-list()
sigmahat2<-matrix(rep(0,N^2),N,N)
for(i in 1:nT){
  Y_e2[[i]]<-Ylist[[i]]-diag(rhohat)%*%W%*%(Ylist[[i]])-Zlist[[i]]%*%betahat
  sigmahat2<-sigmahat2+(Ylist[[i]]-diag(rhohat)%*%W%*%(Ylist[[i]])-Zlist[[i]]%*%betahat)%*%t(Ylist[[i]]-diag(rhohat)%*%W%*%(Ylist[[i]])-Zlist[[i]]%*%betahat)
}
sigmahat2<-sigmahat2/nT
svd(sigmahat2)
Y_e2<-data.frame(matrix(unlist(Y_e2), nrow=nT, byrow=T))

library(psych)
CE4<-corr.test(Y_e2)
CE4$r[upper.tri(CE4$r)]<-0
CE4$p[upper.tri(CE4$p)]<-1

pvalueSigma4<-CE4$p[order(CE4$p)][(N+1):((N+1)*N/2)]
fdr4<-fdrtool(pvalueSigma4,statistic="pvalue")
sum(round(fdr4$qval,5)<0.01)
sum(round(fdr4$qval,5)<0.05)

