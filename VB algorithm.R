rm(list=ls(all=TRUE))
set.seed(4)

setwd('U:\\independent studies\\variational bayes\\binomial LDA') 
source('VB functions.R')
tmp=data.matrix(read.csv('fake data 7.csv',as.is=T))
dat=aggregate.data(tmp)

nobs=5
nloc=nrow(dat)
nspp=ncol(dat)
ncommun=10

#initial values
m0=m1=array(1/ncommun,dim=c(nloc,nspp,ncommun),
            dimnames=list(paste('loc',1:nloc,sep=''),
                          paste('spp',1:nspp,sep=''),
                          paste('comm',1:ncommun,sep='')))
a=b=matrix(1,nloc,ncommun)
c=d=matrix(1,ncommun,nspp)

#priors
gamma=0.1
a.phi=0.01; b.phi=0.99

max.iter=1000
elbo=rep(NA,max.iter)
i=1
thresh=0.0001
delta.elbo=Inf

while (i < max.iter & delta.elbo>thresh){
  print(i)
  tmp=get.ab(ncommun,nloc,nspp,nobs,gamma,m1,m0,dat)
  a=tmp$a
  b=tmp$b
  
  tmp=get.m1m0(nloc,ncommun,nspp,a,b,c,d)
  m1=tmp$m1
  m0=tmp$m0
  
  tmp=get.cd(ncommun,nspp,a.phi,b.phi,nobs,dat,m1,m0)
  c=tmp$c
  d=tmp$d
  
  elbo[i]=get.elbo(a,b,c,d,m1,m0,dat,nobs,nloc,nspp,ncommun,gamma,a.phi,b.phi)
  if (i!=1) delta.elbo=elbo[i]-elbo[i-1]
  i=i+1
}

plot(elbo)

theta=get.theta(a,b)
plot(NA,NA,type='l',ylim=c(0,1),xlim=c(0,nloc))
for (i in 1:ncommun) lines(theta[,i],col=i)

boxplot(theta)