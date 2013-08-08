
NPMLEcmprsk<-function(DATA,C,iteration){

DATA[which(DATA[,2]==0),2]=rep(max(DATA[,2])+1,length(which(DATA[,2]==0)))

DATA=DATA[do.call(order,as.list(as.data.frame(DATA))),]

zdim=dim(DATA)[2]-2
kdim=length(table(DATA[,2]))-1

Data.X=matrix(DATA[,1])
Data.D=matrix(DATA[,2])
Data.Z=as.matrix(DATA[,3:(2+zdim)])

######################## initial values ####################### 

N=length(DATA[,1])
basis=cbind(matrix(1,N,1),Data.Z)

Dk=sapply(1:(kdim+1),function(j) {M.temp<-rep(1,N)
                                  temp=which(Data.D!=j)
                                  M.temp[temp]=rep(0,length(temp))
                                  return(M.temp)})

coe=matrix(0,zdim+1,kdim-1)

rownames(coe)=c("(Intercept)",1:zdim)

b=matrix(0,zdim,kdim)
L=matrix(Data.X,N,kdim)

pLikelihood=function(theta,L){
                              plike.coe=matrix(theta[1:(zdim+1)*(kdim-1)],zdim+1,kdim-1)
                              plike.b=matrix(theta[((zdim+1)*(kdim-1)+1):((zdim+1)*(kdim-1)+zdim*kdim)],zdim,kdim)
                              plike.L=L
                              
                              temp1=exp(basis %*% plike.coe)
                              temp2=temp1/(1+rowSums(temp1))
                              plike.alp=cbind(temp2,1-rowSums(temp2))
                               
                              for (loop in 1:10){
                                                 temp1=plike.alp*exp(-exp(Data.Z %*% plike.b)*plike.L)
                                                 V=temp1/rowSums(temp1)
                                                 
                                                 W=apply((Dk[,1:kdim]+Dk[,kdim+1]*V)*exp(Data.Z %*% plike.b),2,function(list) rev(cumsum(rev(list))))
                                                 W=sapply(1:kdim,function(i) {
                                                                              temp=which(W[,i]==0)
                                                                              if(length(temp)!=0)
                                                                              W[temp,i]=100
                                                                              return(W[,i])})
                                                 
                                                 plike.L=apply(Dk[,1:kdim]/W,2,cumsum)
                                                 }
                              
                              lambda<-plike.L-rbind(rep(0,kdim),plike.L[-N,])
                              temp1=plike.alp*exp(-exp(Data.Z %*% plike.b)*plike.L)
                              temp2=temp1*lambda*exp(Data.Z %*% plike.b)
                              
                              sum(log(temp2^Dk[,1:kdim]))+sum(log(rowSums(temp1)^Dk[,kdim+1]))
                              }

for (i in 1:iteration){

######################## estimate coe ###########################  
  
  temp1=exp(basis %*% coe)
  temp2=temp1/(1+rowSums(temp1))
  alp=cbind(temp2,1-rowSums(temp2))
  
  temp1=alp*exp(-exp(Data.Z %*% b)*L)
  V=temp1/rowSums(temp1)

  temp1=t(basis) %*% (alp[,kdim]*(Dk[,1:(kdim-1)]+Dk[,kdim+1]*V[,1:(kdim-1)]))+C
  temp2=t(basis) %*% (alp[,1:(kdim-1)]*(Dk[,kdim]+Dk[,kdim+1]*V[,kdim]))+C
  coe=coe+log(temp1/temp2)
  
########################estimate beta##########################

  temp1=exp(basis %*% coe)
  temp2=temp1/(1+rowSums(temp1))
  alp=cbind(temp2,1-rowSums(temp2))
  
  temp1=alp*exp(-exp(Data.Z %*% b)*L)
  V=temp1/rowSums(temp1)

  temp1=t(Data.Z) %*% Dk[,1:kdim]
  temp2=t(Data.Z) %*% (L*exp(Data.Z %*% b)*(Dk[,1:kdim]+Dk[,kdim+1]*V[,1:kdim]))
  b=b+log(temp1/temp2)
  
########################estimate lambda#########################
  
  temp1=alp*exp(-exp(Data.Z %*% b)*L)
  V=temp1/rowSums(temp1)

  W=apply((Dk[,1:kdim]+Dk[,kdim+1]*V)*exp(Data.Z%*%b),2,function(list) rev(cumsum(rev(list))))
  
  W=sapply(1:kdim,function(i) {
                               temp=which(W[,i]==0)
                               if(length(temp)!=0)
                               W[temp,i]=100
                               return(W[,i])})
  L=apply(Dk[,1:kdim]/W,2,cumsum)
  }

Sigma.entryij=function(i,j) pLikelihood(c(coe,b)+diag(length(c(coe,b)))[i,]/sqrt(N)+diag(length(c(coe,b)))[,j]/sqrt(N),L)
Sigma.entryi=function(i) pLikelihood(c(coe,b)+diag(length(c(coe,b)))[i,]/sqrt(N),L)

V.Sigma.entryij=Vectorize(Sigma.entryij)

sigmaij=outer(seq(c(coe,b)),seq(c(coe,b)),V.Sigma.entryij)
sigmai=sapply(seq(c(coe,b)),Sigma.entryi)

sigma=-sigmaij+matrix(sigmai,length(c(coe,b)),length(c(coe,b)))+t(matrix(sigmai,length(c(coe,b)),length(c(coe,b))))-matrix(pLikelihood(c(coe,b),L),length(c(coe,b)),length(c(coe,b)))

SD<-sqrt(diag(solve(sigma)))/sqrt(N)
                              
list(Lambda=L,alpha=coe,"SD.of.alpha"=matrix(SD[1:(zdim+1)*(kdim-1)],zdim+1,kdim-1),beta=b,"SD.of.beta"=matrix(SD[((zdim+1)*(kdim-1)+1):((zdim+1)*(kdim-1)+zdim*kdim)],zdim,kdim))

}

