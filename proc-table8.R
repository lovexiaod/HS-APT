rm(list=ls(all=TRUE))
library(MASS);library(Matrix)
library(xtable)
library(plyr)
library(rgenoud)
library(fdrtool)
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


#HS-APT
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

library(xtable)

s<-rep(0,N)
for (l in 1:N){
  s[l]<-(b3[[l]]-2*t(rhohat[l])%*%b1[[l]]-2*t(b22[[l]])%*%b2[[l]]+t(rhohat[l])%*%A11[[l]]%*%rhohat[l]+2*t(b22[[l]])%*%A21[[l]]%*%rhohat[l]+t(b22[[l]])%*%A22[[l]]%*%b22[[l]])/(nT)
}


###significant test
Y_e2<-list()
sigmahat2<-matrix(rep(0,N^2),N,N)
for(i in 1:nT){
  Y_e2[[i]]<-Ylist[[i]]-diag(rhohat)%*%W%*%(Ylist[[i]])-Zlist[[i]]%*%betahat
  sigmahat2<-sigmahat2+(Ylist[[i]]-diag(rhohat)%*%W%*%(Ylist[[i]])-Zlist[[i]]%*%betahat)%*%t(Ylist[[i]]-diag(rhohat)%*%W%*%(Ylist[[i]])-Zlist[[i]]%*%betahat)
}
sigmahat2<-sigmahat2/nT
Y_e2<-data.frame(matrix(unlist(Y_e2), nrow=nT, byrow=T))

library(psych)
CE2<-corr.test(Y_e2)
CE2$r[upper.tri(CE2$r)]<-0
CE2$p[upper.tri(CE2$p)]<-1
round(CE2$r,4)
round(CE2$p,4)
cormatrix<-round(CE2$r,4)
pmatrix<-round(CE2$p,4)
rownames(pmatrix)<-c("Austria","Belgium","Finland","France","Germany","Greece","Ireland","Italy","Netherlands","Portugal","Spain")
colnames(pmatrix)<-c("Austria","Belgium","Finland","France","Germany","Greece","Ireland","Italy","Netherlands","Portugal","Spain")
rownames(cormatrix)<-c("Austria","Belgium","Finland","France","Germany","Greece","Ireland","Italy","Netherlands","Portugal","Spain")
colnames(cormatrix)<-c("Austria","Belgium","Finland","France","Germany","Greece","Ireland","Italy","Netherlands","Portugal","Spain")


Ymatrix<-data.frame(matrix(unlist(Ylist), nrow=nT, byrow=T))
CE<-corr.test(Ymatrix)
CE$r[upper.tri(CE$r)]<-0
CE$p[upper.tri(CE$p)]<-1



round(CE$r,4)
xtable(cormatrix,digits = 4)

