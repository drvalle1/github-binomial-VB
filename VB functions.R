aggregate.data=function(dat){
  ind=which(colnames(dat)=='loc.id')
  locid=dat[,ind]
  dat1=dat[,-ind]
  
  nloc=max(locid)
  nspp=ncol(dat1)
  res=matrix(NA,nloc,nspp)
  for (i in 1:nloc){
    cond=locid==i
    res[i,]=colSums(dat1[cond,])
  }
  res
}
#-----------------------------------
get.ab=function(ncommun,nloc,nspp,nobs,gamma,m1,m0,dat){
  a=b=matrix(NA,nloc,ncommun)
  
  for (i in 1:ncommun){
    nome=paste('comm',i,sep='')
    tmp=dat*m1[,,nome]+(nobs-dat)*m0[,,nome]
    a[,i]=rowSums(tmp)+1
    
    #get cumsum of m1 and m0
    if (i==ncommun) {
      nome=paste('comm',i,sep='')
      cs.m1=m1[,,nome]
      cs.m0=m0[,,nome]
    }
    if (i!=ncommun){
      cs.m1=cs.m0=matrix(0,nloc,nspp)
      for (j in (i+1):ncommun){
        nome=paste('comm',j,sep='')
        cs.m1=cs.m1+m1[,,nome]
        cs.m0=cs.m0+m0[,,nome]
      }  
    }
    tmp=dat*cs.m1+(nobs-dat)*cs.m0
    b[,i]=rowSums(tmp)+gamma
  }
  list(a=a,b=b)
}
#-----------------------------------
get.cd=function(ncommun,nspp,a.phi,b.phi,nobs,dat,m1,m0){
  c=d=matrix(NA,ncommun,nspp)
  for (i in 1:ncommun){
    nome=paste('comm',i,sep='')
    tmp=dat*m1[,,nome]
    c[i,]=colSums(tmp)+a.phi
    
    tmp=(nobs-dat)*m0[,,nome]
    d[i,]=colSums(tmp)+b.phi
  }
  list(c=c,d=d)
}
#----------------------------------
get.m1m0=function(nloc,ncommun,nspp,a,b,c,d){
  #get stb prior piece
  stb=digamma(a)-digamma(a+b)
  stb[,ncommun]=0 #because V_{lK}=1
  
  res=rep(0,nloc)
  for (i in 2:ncommun){
    res=res+digamma(b[,i-1])-digamma(a[,i-1]+b[,i-1])
    stb[,i]=stb[,i]+res
  }

  #get data part
  pdat1=digamma(c)-digamma(c+d)
  pdat0=digamma(d)-digamma(c+d)
  
  m0=m1=array(NA,dim=c(nloc,nspp,ncommun),
              dimnames=list(paste('loc',1:nloc,sep=''),
                            paste('spp',1:nspp,sep=''),
                            paste('comm',1:ncommun,sep='')))
  for (i in 1:nloc){
    for (j in 1:nspp){
      m1.tmp=exp(pdat1[,j]+stb[i,])
      m1[i,j,]=m1.tmp/sum(m1.tmp)
      
      m0.tmp=exp(pdat0[,j]+stb[i,])
      m0[i,j,]=m0.tmp/sum(m0.tmp)
    }
  }
  list(m1=m1,m0=m0)
}
#------------------------------------
get.theta=function(a,b){
  theta=matrix(NA,nloc,ncommun)
  medias=a/(a+b)
  medias[,ncommun]=1
  theta[,1]=medias[,1]
  theta[,2]=medias[,2]*(1-medias[,1])
  for (i in 3:ncommun){
    theta[,i]=medias[,i]*apply(1-medias[,1:(i-1)],1,prod)
  }
  theta
}
#----------------------------------
get.elbo=function(a,b,c,d,m1,m0,dat,nobs,nloc,nspp,ncommun,gamma,a.phi,b.phi){
  #uteis
  dig.a=digamma(a)
  dig.b=digamma(b)
  dig.c=digamma(c)
  dig.d=digamma(d)
  dig.soma.ab=digamma(a+b)
  dig.soma.cd=digamma(c+d)
  prop=dat/nobs
  absent=nobs-dat
  
  #Part 1
  res1=res0=matrix(NA,nloc,nspp)
  for (i in 1:nloc){
    for (j in 1:nspp){
      res1[i,j]=sum(m1[i,j,]*(dig.c[,j]-dig.soma.cd[,j]))
      res0[i,j]=sum(m0[i,j,]*(dig.d[,j]-dig.soma.cd[,j]))
    }
  }
  p1=sum(dat*res1+absent*res0)

  #Part 2

  #get digamma calcs
  pos=dig.a-dig.soma.ab
  neg=dig.b-dig.soma.ab
  
  res=matrix(NA,nloc,ncommun)
  res[,1]=pos[,1]
  res[,2]=pos[,2]+neg[,1]
  for (i in 3:ncommun){ #stick-breaking structure
    res[,i]=pos[,i]+rowSums(neg[,1:(i-1)])
  }
  
  p2=0  
  for (i in 1:ncommun){
    p21=m1[,,i]*dat+m0[,,i]*absent
    p22=matrix(res[,i],nloc,nspp)
    p2=p2+sum(p21*p22)
  }
  
  #Part 3
  tmp=(gamma-1)*(dig.b-dig.soma.ab)
  p3=sum(tmp)
  
  #Part 4
  p41=(a.phi-1)*sum(dig.c-dig.soma.cd)
  p42=(b.phi-1)*sum(dig.d-dig.soma.cd)
  p4=p41+p42
  
  #Part 5
  p5=0
  for (i in 1:ncommun){
    p51=dat*m1[,,i]*log(m1[,,i])  
    p52=absent*m0[,,i]*log(m0[,,i])
    p5=p5+sum(p51)+sum(p52)
  }
  
  #Part 6
  p61=lgamma(a+b)-lgamma(a)-lgamma(b)
  p62=(a-1)*(dig.a-dig.soma.ab)+(b-1)*(dig.b-dig.soma.ab)
  p6=sum(p61)+sum(p62)
  
  #Part 7
  p71=lgamma(c+d)-lgamma(c)-lgamma(d)
  p72=(c-1)*(dig.c-dig.soma.cd)+(d-1)*(dig.d-dig.soma.cd)
  p7=sum(p71)+sum(p72)
  
  p1+p2+p3+p4-p5-p6-p7
}