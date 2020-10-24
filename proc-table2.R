rm(list=ls(all=TRUE))
library(MASS);library(Matrix)
library(xtable)
library(plyr)
library(rgenoud)
GI<-read.csv("G:/joint manuscript/GI.csv")
colnames(GI)[1] <- 'Month'
Mom<-read.csv("G:/joint manuscript/Mom.csv")
GI<-cbind(GI,rep(Mom[,2],each=11))
tdata<-dlply(GI,.(Month))#
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
c1<-c(1:5);c2<-c(1:3);c3<-c(1:3,6);c4<-c(1:4);c5<-c(1:4,6);c6<-c(1:6) ## c1:c6 denotes 6 different combination
for (i in 1:nT){
  Xlist[[i]]=Xlist[[i]][,c2]
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


library(xtable)
xtable(W,digits=3)


#S-APT
Wt<-list();Wt2<-list()
for (l in 1:N){
  Wt[[l]]<-list()
}
for (l in 1:N){
  Wt2[[l]]<-list()
}
A11<-list();A12<-list();A21<-list();A22<-list();b1<-list();b2<-list();b3<-list()
A112<-list();A122<-list();A212<-list();A222<-list();b12<-list();b22<-list();b32<-list()

for (l in 1:N){
  A11[[l]]<-0;A12[[l]]<-0;A21[[l]]<-0;A22[[l]]<-0;b1[[l]]<-0;b2[[l]]<-0;b3[[l]]<-0
  for (i in 1:nT){
    Wt[[l]][[i]]<-W[l,]*Ylist[[i]] #Nx1
    A11[[l]]<-A11[[l]]+Wt[[l]][[i]]%*%t(Wt[[l]][[i]]) #NxN
    A21[[l]]<-A21[[l]]+Xlist[[i]][l,]%*%t(Wt[[l]][[i]]) #PxN
    A22[[l]]<-A22[[l]]+Xlist[[i]][l,]%*%t(Xlist[[i]][l,]) #PxP
    b1[[l]]<-b1[[l]]+Wt[[l]][[i]]*Ylist[[i]][l] #Nx1
    b2[[l]]<-b2[[l]]+Xlist[[i]][l,]*Ylist[[i]][l] #Px1
    b3[[l]]<-b3[[l]]+Ylist[[i]][[l]]*Ylist[[i]][[l]] #1x1
  }
  A112[[l]]<-0;A122[[l]]<-0;A212[[l]]<-0;A222[[l]]<-0;b12[[l]]<-0;b22[[l]]<-0;b32[[l]]<-0
  for (i in 1:nT){
    Wt2[[l]][[i]]<-W[l,]*Ylist[[i]] #Nx1
    A112[[l]]<-A112[[l]]+Wt2[[l]][[i]]%*%t(Wt2[[l]][[i]]) #NxN
    A212[[l]]<-A212[[l]]+X1list[[i]][l,]%*%t(Wt2[[l]][[i]]) #PxN
    A222[[l]]<-A222[[l]]+X1list[[i]][l,]%*%t(X1list[[i]][l,]) #PxP
    b12[[l]]<-b12[[l]]+Wt2[[l]][[i]]*Ylist[[i]][l] #Nx1
    b22[[l]]<-b22[[l]]+X1list[[i]][l,]*Ylist[[i]][l] #Px1
    b32[[l]]<-b32[[l]]+Ylist[[i]][[l]]*Ylist[[i]][[l]] #1x1
  }
}

DH_H<-function(rho){
  b<-list()
  for (l in 1:N){
    b[[l]]<-solve(A22[[l]])%*%(b2[[l]]-A21[[l]]%*%rep(rho,N)) #Px1
  }
  s<-rep(0,N)
  for (l in 1:N){
    s[l]<-(b3[[l]]-2*t(rep(rho,N))%*%b1[[l]]-2*t(b[[l]])%*%b2[[l]]+t(rep(rho,N))%*%A11[[l]]%*%rep(rho,N)+2*t(b[[l]])%*%A21[[l]]%*%rep(rho,N)+t(b[[l]])%*%A22[[l]]%*%b[[l]])/(nT)
  }
  lc<-(-N*nT*log(2*pi))/2-nT*sum(log(det(diag(s))))/2+nT*log(det((diag(1,N)-t(diag(rep(rho,N)))%*%t(W))%*%(diag(1,N)-W%*%diag(rep(rho,N)))))/2-N*nT/2
  #  lc<--(N*(nT))*log(2*pi*s)/2+(nT)*log(det((diag(1,N)-rho*t(W))%*%(diag(1,N)-W*rho)))/2-N*(nT)/2
  return(lc)
}
parlist<-list()
betalist<-list()
likelihoodlist<-list()
threepartsofL<-list()
nrep2<-10
set.seed(1234)
for (i in 1:nrep2){
  print(c(i,i))
  g_H<-genoud(DH_H, nvars=1, max=TRUE, pop.size=1000,print.level = F)
  likelihoodlist[[i]]<-g_H$value
  parlist[[i]]<-g_H$par
  likelihoodlist[[i]]<-g_H$value
}
rho_ls<-parlist[[which.max(likelihoodlist)]]

b<-list()
for (l in 1:N){
  b[[l]]<-solve(A22[[l]])%*%(b2[[l]]-A21[[l]]%*%rep(rho_ls,N)) 
}

beta_ls<-rep(0,N*p)
for (l in 1:N){
  beta_ls[(1+(l-1)*p):(l*p)]<-b[[l]]
}

s<-rep(0,N)
for (l in 1:N){
  s[l]<-(b3[[l]]-2*t(rep(rho_ls,N))%*%b1[[l]]-2*t(b[[l]])%*%b2[[l]]+t(rep(rho_ls,N))%*%A11[[l]]%*%rep(rho_ls,N)+2*t(b[[l]])%*%A21[[l]]%*%rep(rho_ls,N)+t(b[[l]])%*%A22[[l]]%*%b[[l]])/(nT)
}

AIC_MS_11<--2*g_H$value+2*(N*(p+1)+1)
BIC_MS_11<--2*g_H$value+log(N*(nT))*(N*(p+1)+1)
AIC_MS_11
BIC_MS_11
L1<-likelihoodlist[[which.max(likelihoodlist)]]


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
set.seed(1234)
for (i in 1:nrep2){
  print(c(i,i))
  g_H<-genoud(DH_H2, nvars=N, max=TRUE, pop.size=1000,print.level = F)
  likelihoodlist2[[i]]<-g_H$value
  parlist2[[i]]<-g_H$par
  likelihoodlist2[[i]]<-g_H$value
}

rho_ls<-parlist2[[which.max(likelihoodlist2)]]

b<-list()
for (l in 1:N){
  b[[l]]<-solve(A22[[l]])%*%(b2[[l]]-A21[[l]]%*%rho_ls[l]) 
}

beta_ls<-rep(0,N*p)
for (l in 1:N){
  beta_ls[(1+(l-1)*p):(l*p)]<-b[[l]]
}

s<-rep(0,N)
for (l in 1:N){
  s[l]<-(b3[[l]]-2*t(rho_ls[l])%*%b1[[l]]-2*t(b[[l]])%*%b2[[l]]+t(rho_ls[l])%*%A11[[l]]%*%rho_ls[l]+2*t(b[[l]])%*%A21[[l]]%*%rho_ls[l]+t(b[[l]])%*%A22[[l]]%*%b[[l]])/(nT)
}

AIC_HMS_11list<-list()
for (i in 1:nrep2){
  AIC_HMS_11list[[i]]<--2*likelihoodlist2[[i]]+2*(N*(p+2))
}
AIC_HMS_11<--2*likelihoodlist2[[which.max(likelihoodlist2)]]+2*(N*(p+2))
BIC_HMS_11<--2*likelihoodlist2[[which.max(likelihoodlist2)]]+log(N*(nT))*(N*(p+1)+1)
AIC_HMS_11
BIC_HMS_11

L2<-likelihoodlist2[[which.max(likelihoodlist2)]]
2*L2-2*L1
1-pchisq(2*L2-2*L1,N-1)

