rm(list=ls(all=TRUE))
library(MASS);library(Matrix)
library(xtable)
library(plyr)
library(rgenoud)
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



Wt<-list()
for (l in 1:N){
  Wt[[l]]<-list()
}
Wt_1<-list()
for (l in 1:N){
  Wt_1[[l]]<-list()
}
A11<-list();A12<-list();A21<-list();A22<-list();b1<-list();b2<-list();b3<-list()
A11_1<-list();A12_1<-list();A21_1<-list();A22_1<-list();b1_1<-list();b2_1<-list();b3_1<-list()
for (l in 1:N){
  A11[[l]]<-0;A12[[l]]<-0;A21[[l]]<-0;A22[[l]]<-0;b1[[l]]<-0;b2[[l]]<-0;b3[[l]]<-0
  A11_1[[l]]<-0;A12_1[[l]]<-0;A21_1[[l]]<-0;A22_1[[l]]<-0;b1_1[[l]]<-0;b2_1[[l]]<-0;b3_1[[l]]<-0
  for (i in 1:nT){
    Wt[[l]][[i]]<-W[l,]%*%Ylist[[i]] #1x1
    A11[[l]]<-A11[[l]]+(Wt[[l]][[i]])^2 #1x1
    A21[[l]]<-A21[[l]]+Xlist[[i]][l,]%*%Wt[[l]][[i]] #Px1
    A22[[l]]<-A22[[l]]+Xlist[[i]][l,]%*%t(Xlist[[i]][l,]) #PxP same
    b1[[l]]<-b1[[l]]+Ylist[[i]][l]*Wt[[l]][[i]] #1x1
    b2[[l]]<-b2[[l]]+Xlist[[i]][l,]*Ylist[[i]][l] #Px1 same
    b3[[l]]<-b3[[l]]+Ylist[[i]][[l]]*Ylist[[i]][[l]] #1x1 same
    Wt_1[[l]][[i]]<-W[l,]%*%Ylist[[i]] #1x1
    A11_1[[l]]<-A11_1[[l]]+(Wt_1[[l]][[i]])^2 #1x1
    A21_1[[l]]<-A21_1[[l]]+X1list[[i]][l,]%*%Wt_1[[l]][[i]] #Px1
    A22_1[[l]]<-A22_1[[l]]+X1list[[i]][l,]%*%t(X1list[[i]][l,]) #PxP same
    b1_1[[l]]<-b1_1[[l]]+Ylist[[i]][l]*Wt_1[[l]][[i]] #1x1
    b2_1[[l]]<-b2_1[[l]]+X1list[[i]][l,]*Ylist[[i]][l] #Px1 same
    b3_1[[l]]<-b3_1[[l]]+Ylist[[i]][[l]]*Ylist[[i]][[l]] #1x1 same
  }
}


DH_H1<-function(rho){
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

DH_H2<-function(rho_1){
  b_1<-list()
  for (l in 1:N){
    b_1[[l]]<-solve(A22_1[[l]])%*%(b2_1[[l]]-A21_1[[l]]*rho_1[l]) #Px1
  }
  s_1<-rep(0,N)
  for (l in 1:N){
    s_1[l]<-(b3_1[[l]]-2*b1_1[[l]]*rho_1[l]-2*t(b_1[[l]])%*%b2_1[[l]]+(rho_1[l]^2)*A11_1[[l]]+2*t(b_1[[l]])%*%A21_1[[l]]*rho_1[l]+t(b_1[[l]])%*%A22_1[[l]]%*%b_1[[l]])/(nT)
  }
  lc<-(-N*nT*log(2*pi))/2-nT*log(det(diag(s_1)))/2+nT*log(det((diag(1,N)-t(W)%*%t(diag(c(rho_1))))%*%(diag(1,N)-diag(c(rho_1))%*%W)))/2-N*nT/2
  return(lc)
}

parlist<-list()
betalist<-list()
likelihoodlist<-list()
threepartsofL<-list()
parlist_1<-list()
betalist_1<-list()
likelihoodlist_1<-list()
threepartsofL_1<-list()
nrep2<-5

for (i in 1:nrep2){
  print(c(i,i))
  g_H<-genoud(DH_H1, nvars=N, max=TRUE, pop.size=1000,print.level = F)
  g_H_1<-genoud(DH_H2, nvars=N, max=TRUE, pop.size=1000,print.level = F)
  likelihoodlist[[i]]<-g_H$value
  likelihoodlist_1[[i]]<-g_H_1$value
  parlist[[i]]<-g_H$par
  parlist_1[[i]]<-g_H_1$par
}
rho_ls<-parlist[[which.max(likelihoodlist)]]
rho_ls_1<-parlist_1[[which.max(likelihoodlist)]]

b<-list()
b_1<-list()
for (l in 1:N){
  b[[l]]<-solve(A22[[l]])%*%(b2[[l]]-A21[[l]]%*%rho_ls[l]) 
  b_1[[l]]<-solve(A22_1[[l]])%*%(b2_1[[l]]-A21_1[[l]]%*%rho_ls_1[l]) 
}

beta_ls<-rep(0,N*p)
beta_ls_1<-rep(0,N*p)
for (l in 1:N){
  beta_ls[(1+(l-1)*p):(l*p)]<-b[[l]]
  beta_ls_1[(1+(l-1)*(p+1)):(l*(p+1))]<-b_1[[l]]
}

s<-rep(0,N)
s_1<-rep(0,N)
for (l in 1:N){
  s[l]<-(b3[[l]]-2*t(rho_ls[l])%*%b1[[l]]-2*t(b[[l]])%*%b2[[l]]+t(rho_ls[l])%*%A11[[l]]%*%rho_ls[l]+2*t(b[[l]])%*%A21[[l]]%*%rho_ls[l]+t(b[[l]])%*%A22[[l]]%*%b[[l]])/(nT)
  s_1[l]<-(b3_1[[l]]-2*t(rho_ls_1[l])%*%b1_1[[l]]-2*t(b_1[[l]])%*%b2_1[[l]]+t(rho_ls_1[l])%*%A11_1[[l]]%*%rho_ls_1[l]+2*t(b_1[[l]])%*%A21_1[[l]]%*%rho_ls_1[l]+t(b_1[[l]])%*%A22_1[[l]]%*%b_1[[l]])/(nT)
}


rho_ls
beta_ls
s
