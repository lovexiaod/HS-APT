rm(list=ls(all=TRUE))#清除R工作环境中的全部东西
library(MASS);library(psych)
library(xtable);library(Matrix)
library(plyr)
library(rgenoud)

GI<-read.csv("G:/joint manuscript/GI.csv")
colnames(GI)[1] <- 'Month'

tdata<-dlply(GI,.(Month))#dataframe转list
tdata<-tdata[(length(tdata)-137-24):length(tdata)]
#tdata<-tdata[-(13:24)]

#1MKT,2SMB,3HML,4RMW,5CMA1,3,5,6
cc<-list()
cc[[1]]<-c(1:5);cc[[2]]<-c(1,4,5,6);cc[[3]]<-c(1,2,3)
cc[[4]]<-c(1,3,4,5,6);cc[[5]]<-c(1,2,4,5,6);cc[[6]]<-c(1,2,3,4,6)
cc[[7]]<-c(1,2,3,5,6);cc[[8]]<-c(1,3,4,6);cc[[9]]<-c(1,2,4,6)
cc[[10]]<-c(1,3,5,6);cc[[11]]<-c(1,2,5,6);cc[[12]]<-c(1:6);cc[[13]]<-c(1:3,6)
RE<-matrix(rep(0,8*length(cc)),8,length(cc))
rhoslm<-list();betaslm<-list()

Xlist<-list()
Ylist<-list()

nT=length(tdata)
N=11

for (i in 1:nT){
  #Ylist[[i]]=tdata[[i]][,3]
  Xlist[[i]]=as.matrix(cbind(tdata[[i]][,c(4:8)]))
  Ylist[[i]]=tdata[[i]][,3]-tdata[[i]][,9]
}

for (i in 1:nT){
  Xlist[[i]]=Xlist[[i]][,cc[[1]]]
}
p=ncol(Xlist[[1]])

XBlist<-list()
for (i in 1:nT){
  XBlist[[i]]=cbind(diag(1,N),Xlist[[i]])
}
Zlist<-list()
for (i in 1:nT){
  xx<-list()
  for (j in 1:11){
    xx[[j]]<-t(c(Xlist[[i]][j,]))
  }
  Zlist[[i]]<-as.matrix(bdiag(xx))
}
Z1list<-list()
for (i in 1:nT){
  xx<-list()
  for (j in 1:11){
    xx[[j]]<-t(c(1,Xlist[[i]][j,]))
  }
  Z1list[[i]]<-as.matrix(bdiag(xx))
}
for (i in 1:nT){
  Ylist[[i]]=tdata[[i]][,3]-tdata[[i]][,9]
}
sumY<-rep(0,N)
for (i in 1:nT){
  sumY<-sumY+Ylist[[i]]
}
YC<-sumY/nT

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

#APT
Wt1<-list();BA11<-0;BA12<-0;BA21<-0;BA22<-0;Bb1<-0;Bb2<-0;Bb3<-0
Wt1_1<-list();BA11_1<-0;BA12_1<-0;BA21_1<-0;BA22_1<-0;Bb1_1<-0;Bb2_1<-0;Bb3_1<-0
for (i in 1:nT){
  Wt1[[i]]<-W%*%Ylist[[i]]
  BA11<-BA11+t(Wt1[[i]])%*%Wt1[[i]]
  BA12<-BA12+t(Wt1[[i]])%*%Zlist[[i]]
  BA21<-BA21+t(Zlist[[i]])%*%Wt1[[i]]
  BA22<-BA22+t(Zlist[[i]])%*%Zlist[[i]]
  Bb1<-Bb1+t(Wt1[[i]])%*%Ylist[[i]]
  Bb2<-Bb2+t(Zlist[[i]])%*%Ylist[[i]]
  Bb3<-Bb3+t(Ylist[[i]])%*%Ylist[[i]]
  Wt1_1[[i]]<-W%*%Ylist[[i]]
  BA11_1<-BA11_1+t(Wt1_1[[i]])%*%Wt1_1[[i]]
  BA12_1<-BA12_1+t(Wt1_1[[i]])%*%Z1list[[i]]
  BA21_1<-BA21_1+t(Z1list[[i]])%*%Wt1_1[[i]]
  BA22_1<-BA22_1+t(Z1list[[i]])%*%Z1list[[i]]
  Bb1_1<-Bb1_1+t(Wt1_1[[i]])%*%Ylist[[i]]
  Bb2_1<-Bb2_1+t(Z1list[[i]])%*%Ylist[[i]]
  Bb3_1<-Bb3_1+t(Ylist[[i]])%*%Ylist[[i]]
}
LL_S<-function(rho){
  b<-solve(BA22)%*%(Bb2-BA21*rho)
  s<-(Bb3-2*rho*Bb1-2*t(b)%*%Bb2+rho*BA11*rho+2*t(b)%*%BA21*rho+t(b)%*%BA22%*%b)/(N*nT)
  lc<--(N*nT)*log(2*pi*s)/2+nT*log(det((diag(1,N)-rho*t(W))%*%(diag(1,N)-W*rho)))/2-N*nT/2
  return(lc)
}
LL_S_2<-function(rho){
  b<-solve(BA22_1)%*%(Bb2_1-BA21_1*rho)
  s<-(Bb3_1-2*rho*Bb1_1-2*t(b)%*%Bb2_1+rho*BA11_1*rho+2*t(b)%*%BA21_1*rho+t(b)%*%BA22_1%*%b)/(N*nT)
  lc<--(N*nT)*log(2*pi*s)/2+nT*log(det((diag(1,N)-rho*t(W))%*%(diag(1,N)-W*rho)))/2-N*nT/2
  return(lc)
}

nrep1<-10
g_slist<-list()
likelilist<-rep(0,nrep1)
for (i in 1:nrep1){
  print(i)
  g_slist[[i]]<-genoud(LL_S, nvars=1, max=TRUE, pop.size=1000,print.level = F)
  likelilist[i]<-g_slist[[i]]$value
}
g_S<-g_slist[[which.max(likelilist)]]
g_S_2<-genoud(LL_S_2, nvars=1,max=TRUE, pop.size=1000,print.level = F)
LR_SLM<-1-pchisq(2*abs(g_S$value-g_S_2$value),11)
AIC_SLM_1<--2*g_S$value+2*(N*p+1)
AIC_SLM_2<--2*g_S_2$value+2*(N+N*p+1)
BIC_SLM_1<--2*g_S$value+log(N*nT)*(N*p+1)
BIC_SLM_2<--2*g_S_2$value+log(N*nT)*(N+N*p+1)


rho_slm<-g_S$par;beta_slm<-solve(BA22)%*%(Bb2-BA21%*%g_S$par)
rho_slm_2<-g_S_2$par;beta_slm_2<-solve(BA22_1)%*%(Bb2_1-BA21_1%*%g_S_2$par)
rhoslm<-rho_slm;betaslm<-beta_slm
SSE_slm<-0
SST_slm<-0
SSE_slm_2<-0
for (i in 1:nT){
  SSE_slm<-SSE_slm+((diag(1,N)-as.numeric(rho_slm)*W)%*%Ylist[[i]]-Zlist[[i]]%*%beta_slm)^2
  SST_slm<-SST_slm+(Ylist[[i]]-YC)^2
  SSE_slm_2<-SSE_slm_2+(Ylist[[i]]-as.numeric(rho_slm_2)*W%*%Ylist[[i]]-(Z1list[[i]]%*%beta_slm_2))^2
}
R2_SLM<-1-SSE_slm/SST_slm
R2_SLM_2<-1-SSE_slm_2/SST_slm
ADR2_SLM<-sum(1-(1-R2_SLM)*(nT-1)/(nT-N*p-1))/N
ADR2_SLM_2<-sum(1-(1-R2_SLM_2)*(nT-1)/(nT-N*p-N-1))/N

Errlist<-list()
Errmatrix<-matrix(rep(0,nT*N),N,nT)
for (i in 1:nT){
  #kronecker(diag(rep(1,N)),t(Xlist[[i]][1,]))
  Errlist[[i]]<-Ylist[[i]]-rho_slm*W%*%Ylist[[i]]-Zlist[[i]]%*%beta_slm
  Errmatrix[,i]<-Errlist[[i]]
}

PCAerrtest1<-princomp(t(Errmatrix)) 
PCAerr<-princomp(t(Errmatrix), cor = TRUE)
PCA1<-loadings(PCAerrtest1)[,1]


###SAR model with first PCA vector about error term
EXlist<-list()
for (i in 1:nT){
  EXlist[[i]]<-PCA1*(sum(PCA1*Errlist[[i]]))
}


Wt1<-list();BA11<-0;BA12<-0;BA21<-0;BA22<-0;Bb1<-0;Bb2<-0;Bb3<-0
Wt1_1<-list();BA11_1<-0;BA12_1<-0;BA21_1<-0;BA22_1<-0;Bb1_1<-0;Bb2_1<-0;Bb3_1<-0
ZElist<-list()
for (i in 1:nT){
  xx<-list()
  for (j in 1:11){
    xx[[j]]<-t(c(EXlist[[i]][j]))
  }
  ZElist[[i]]<-as.matrix(bdiag(xx))
}
Z1Elist<-list()
for (i in 1:nT){
  xx<-list()
  for (j in 1:11){
    #xx[[j]]<-t(c(1,EXlist[[i]][j,]))
    xx[[j]]<-t(c(1,EXlist[[i]][j]))
  }
  Z1Elist[[i]]<-as.matrix(bdiag(xx))
}

for (i in 1:nT){
  Wt1[[i]]<-W%*%Errlist[[i]]
  BA11<-BA11+t(Wt1[[i]])%*%Wt1[[i]]
  BA12<-BA12+t(Wt1[[i]])%*%ZElist[[i]]
  BA21<-BA21+t(ZElist[[i]])%*%Wt1[[i]]
  BA22<-BA22+t(ZElist[[i]])%*%ZElist[[i]]
  Bb1<-Bb1+t(Wt1[[i]])%*%Errlist[[i]]
  Bb2<-Bb2+t(ZElist[[i]])%*%Errlist[[i]]
  Bb3<-Bb3+t(Errlist[[i]])%*%Errlist[[i]]
  Wt1_1[[i]]<-W%*%Errlist[[i]]
  BA11_1<-BA11_1+t(Wt1_1[[i]])%*%Wt1_1[[i]]
  BA12_1<-BA12_1+t(Wt1_1[[i]])%*%Z1Elist[[i]]
  BA21_1<-BA21_1+t(Z1Elist[[i]])%*%Wt1_1[[i]]
  BA22_1<-BA22_1+t(Z1Elist[[i]])%*%Z1Elist[[i]]
  Bb1_1<-Bb1_1+t(Wt1_1[[i]])%*%Errlist[[i]]
  Bb2_1<-Bb2_1+t(Z1Elist[[i]])%*%Errlist[[i]]
  Bb3_1<-Bb3_1+t(Errlist[[i]])%*%Errlist[[i]]
}
LL_SE<-function(rho){
  b<-solve(BA22)%*%(Bb2-BA21*rho)
  s<-(Bb3-2*rho*Bb1-2*t(b)%*%Bb2+rho*BA11*rho+2*t(b)%*%BA21*rho+t(b)%*%BA22%*%b)/(N*nT)
  lc<--(N*nT)*log(2*pi*s)/2+nT*log(det((diag(1,N)-rho*t(W))%*%(diag(1,N)-W*rho)))/2-N*nT/2
  return(lc)
}

LL_S_2E<-function(rho){
  b<-solve(BA22_1)%*%(Bb2_1-BA21_1*rho)
  s<-(Bb3_1-2*rho*Bb1_1-2*t(b)%*%Bb2_1+rho*BA11_1*rho+2*t(b)%*%BA21_1*rho+t(b)%*%BA22_1%*%b)/(N*nT)
  lc<--(N*nT)*log(2*pi*s)/2+nT*log(det((diag(1,N)-rho*t(W))%*%(diag(1,N)-W*rho)))/2-N*nT/2
  return(lc)
}

nrep1<-5
g_slistE<-list()
likelilistE<-rep(0,nrep1)
for (i in 1:nrep1){
  print(i)
  g_slistE[[i]]<-genoud(LL_SE, nvars=1, max=TRUE, pop.size=1000,print.level = F)
  likelilistE[i]<-g_slistE[[i]]$value
}
g_SE<-g_slistE[[which.max(likelilistE)]]
g_S2E<-genoud(LL_S_2E, nvars=1,max=TRUE, pop.size=1000,print.level = F)

Errbetahat1<-solve(BA22)%*%Bb2
Errshat1<-(Bb3-2*t(Errbetahat1)%*%Bb2+t(Errbetahat1)%*%BA22%*%Errbetahat1)/(N*nT)
Errlikelihood1<--(N*nT)*log(2*pi*Errshat1)/2+nT*log(det((diag(1,N))%*%(diag(1,N))))/2-N*nT/2

LR_Err1<-1-pchisq(2*abs(g_SE$value-Errlikelihood1),1)

LR_Err1


###check whether S-APT or HS-APT is more suitable
rm(list=ls(all=TRUE))
library(MASS);library(Matrix)
library(xtable)
library(plyr)
library(rgenoud)
GI<-read.csv("G:/joint manuscript/GI.csv")
colnames(GI)[1] <- 'Month'
Mom<-read.csv("G:/joint manuscript/Mom.csv")
GI<-cbind(GI,rep(Mom[,2],each=11))
tdata<-dlply(GI,.(Month))#dataframe转list
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
  #  lc<--(N*(nT))*log(2*pi*s)/2+(nT)*log(det((diag(1,N)-rho*t(W))%*%(diag(1,N)-W*rho)))/2-N*(nT)/2 同方差
  return(lc)
}
parlist<-list()
betalist<-list()
likelihoodlist<-list()
threepartsofL<-list()
nrep2<-10

for (i in 1:nrep2){
  print(c(i,i))
  g_H<-genoud(DH_H, nvars=1, max=TRUE, pop.size=1000,print.level = F)
  likelihoodlist[[i]]<-g_H$value
  parlist[[i]]<-g_H$par
  #slist[[i]]<-(b3-2*t(parlist[[i]])%*%b1-2*t(betalist[[i]])%*%b2+t(parlist[[i]])%*%A11%*%parlist[[i]]+2*t(betalist[[i]])%*%A21%*%parlist[[i]]+t(betalist[[i]])%*%A22%*%betalist[[i]])/(N*nT)
  #threepartsofL[[i]]<-c(-(N*nT)*log(slist[[i]])/2,nT*log(det((diag(1,N)-t(diag(c(parlist[[i]])))%*%t(W))%*%(diag(1,N)-W%*%diag(c(parlist[[i]])))))/2,-N*nT/2)
  likelihoodlist[[i]]<-g_H$value
}
#rho_ls<-parlist[[which.max(likelihoodlist)]] #to be continued
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
  #threepartsofL[[i]]<-c(-(N*nT)*log(slist[[i]])/2,nT*log(det((diag(1,N)-t(diag(c(parlist[[i]])))%*%t(W))%*%(diag(1,N)-W%*%diag(c(parlist[[i]])))))/2,-N*nT/2)
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

L2<-likelihoodlist2[[which.max(likelihoodlist2)]]
2*L2-2*L1
1-pchisq(2*L2-2*L1,N-1)
