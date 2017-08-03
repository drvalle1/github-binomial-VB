rm(list=ls(all=TRUE))
set.seed(4)

nloc=1000
nobs=5
nspp=90
y=matrix(NA,nobs*nloc,nspp)
loc.id=rep(1:nloc,each=nobs)

ncommun=7

#generate thetas
base=floor(nloc/(ncommun-2))

x=seq(from=-1,to=1,length.out=base)
y=sqrt(1-(x^2))*0.1
min1=0.0001
y[y<min1]=min1
# plot(x,y)

init=floor(nloc/ncommun)
seq1=c(seq(from=1,to=nloc,by=init),nloc)

theta=matrix(min1,nloc,ncommun)
for (i in 1:ncommun){
  seq2=seq1[i]:(seq1[i]+base-1)
  seq3=seq2[seq2<=nloc]
  theta[seq3,i]=y[1:length(seq3)]
}
theta=theta/matrix(apply(theta,1,sum),nloc,ncommun)
theta.true=theta

plot(NA,NA,xlim=c(0,nloc),ylim=c(0,1))
for (i in 1:ncommun) lines(1:nloc,theta[,i],col=i)

#generate phi's
makesparse=matrix(rbinom(ncommun*nspp,size=1,prob=0.3),ncommun,nspp)
phi.true=phi=matrix(runif(ncommun*nspp,min=0.8,max=1),ncommun,nspp)*makesparse

par(mfrow=c(3,2),mar=rep(1,4))
for (i in 1:ncommun) plot(phi[i,],type='h',ylim=c(0,1))

#generate data
probs=theta%*%phi
hist(probs)
probs1=probs[loc.id,]

#number of observations per location
obs=matrix(rbinom(nrow(probs1)*ncol(probs1),size=1,prob=probs1),nrow(probs1),ncol(probs1))
colnames(obs)=paste('spp',1:nspp,sep='')

par(mfrow=c(1,1))
image(obs)

obs1=cbind(obs,loc.id)

setwd('U:\\independent studies\\variational bayes\\binomial LDA') 
nome=paste(c('fake data ','theta '),ncommun,'.csv',sep='')
write.csv(obs1,nome[1],row.names=F)
write.csv(theta,nome[2],row.names=F)
