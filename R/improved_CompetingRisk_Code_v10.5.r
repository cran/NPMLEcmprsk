NPMLEcmprsk<-function(DATA,censoring.coding=0,alpha.stable.parameter=100,beta.stable.parameter=100,initial.alpha=0,initial.beta=0,threshold=0,iteration=5000){

DATA=DATA[do.call(order,as.list(as.data.frame(DATA))),]

if(censoring.coding==0){                                     # re-write by FangYu 20140319
 if ( sum(DATA[,2]==censoring.coding)!=0 ){                  # re-write by FangYu 20140319
  fishmax=max(DATA[,2])+1                                    # re-write by FangYu 20140422
  DATA[which(DATA[,2]==censoring.coding),2]=fishmax          # re-write by FangYu 20140422
  censoring.coding= fishmax }                                # re-write by FangYu 20140422
}                                                            # re-write by FangYu 20140319

zdim=dim(DATA)[2]-2
kname<- names(table(DATA[,2])  )        # re-write by FangYu 20140319
kdim<- sum(kname != censoring.coding)   # re-write by FangYu 20140319

Data.X=matrix(DATA[,1])
Data.D=matrix(DATA[,2])
Data.Z=as.matrix(DATA[,3:(2+zdim)])

# shift<-min(as.matrix(DATA[,3:(2+zdim)]))

# if(shift<0)
# Data.Z=as.matrix(DATA[,3:(2+zdim)])+abs(shift)+0.1


######################## initial values #######################

N=length(DATA[,1])
basis=cbind(matrix(1,N,1),Data.Z)

Dk=sapply(1:(kdim+1),function(j) {M.temp<-rep(1,N)
                                  temp=which(Data.D!=j)
                                  M.temp[temp]=rep(0,length(temp))
                                  return(M.temp)})

coe=matrix(initial.alpha,zdim+1,kdim-1)
b=matrix(initial.beta,zdim,kdim)
# L=matrix(Data.X,N,kdim)
L=matrix(0,N,kdim)

C.coe=matrix(0,zdim+1,kdim-1)
C.b=matrix(0,zdim,kdim)

stop_sum=0

pLikelihood=function(theta,L){
                              plike.L=L
                              
                              plike.coe=matrix(theta[(kdim*zdim+1):(zdim*kdim+(zdim+1)*(kdim-1))],zdim+1,kdim-1)
                              plike.b=matrix(theta[1:(kdim*zdim)],zdim,kdim)
                              
                              temp1=exp(basis %*% plike.coe)
                              temp2=temp1/(1+rowSums(temp1))
                              plike.alp=cbind(temp2,1-rowSums(temp2))

                              for (loop in 1:30){
                                                 temp1=plike.alp*exp(-exp(Data.Z %*% plike.b)*plike.L)
                                                 
                                                 if(sum(rowSums(temp1)==0)!=0){
                                                                               message("The value of W is zero, numerical approximation is applied.")  
                                                                               ## using log transformation
                                                                               temp1=log(plike.alp)-exp(Data.Z %*% plike.b)*plike.L
                                                                               }
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

Likelihood=function(theta,L){
                             plike.L=L
                              
                             plike.coe=matrix(theta[(kdim*zdim+1):(zdim*kdim+(zdim+1)*(kdim-1))],zdim+1,kdim-1)
                             plike.b=matrix(theta[1:(kdim*zdim)],zdim,kdim)
                              
                             temp1=exp(basis %*% plike.coe)
                             temp2=temp1/(1+rowSums(temp1))
                             plike.alp=cbind(temp2,1-rowSums(temp2))
                             
                             lambda<-plike.L-rbind(rep(0,kdim),plike.L[-N,])
                             temp1=plike.alp*exp(-exp(Data.Z %*% plike.b)*plike.L)
                             temp2=temp1*lambda*exp(Data.Z %*% plike.b)
                             sum(log(temp2^Dk[,1:kdim]))+sum(log(rowSums(temp1)^Dk[,kdim+1]))
                             }

for (i in 1:iteration){

######################## estimate lambda #########################
 
  temp1.alpha=exp(basis %*% coe)
  temp2.alpha=temp1.alpha/(1+rowSums(temp1.alpha))
  
  alp=cbind(temp2.alpha,1-rowSums(temp2.alpha))
 
  temp=Data.Z %*% b
  temp1=alp*exp(-exp(temp)*L)
    
  if(sum(rowSums(temp1)==0)!=0){
                                message("The value of W is zero, numerical approximation is applied.")  ## using log transformation
                                temp1=log(alp)-exp(temp)*L
                                }
  
  V=temp1/rowSums(temp1)
    
  W=apply((Dk[,1:kdim]+Dk[,kdim+1]*V)*exp(temp),2,function(list) rev(cumsum(rev(list))))
  W=sapply(1:kdim,function(i) {
                               temp=which(W[,i]==0)
                               if(length(temp)!=0)
                               W[temp,i]=100
                               return(W[,i])})
  L=apply(Dk[,1:kdim]/W,2,cumsum)

  
######################## estimate coe ###########################
    
  temp1.alpha=exp(basis %*% coe)
  temp2.alpha=temp1.alpha/(1+rowSums(temp1.alpha))
  
  alp=cbind(temp2.alpha,1-rowSums(temp2.alpha))
  
  temp=Data.Z %*% b
  temp1.alpha=alp*exp(-exp(temp)*L)
    
  if(sum(rowSums(temp1.alpha)==0)!=0){
                                message("The value of W is zero, numerical approximation is applied.")  ## using log transformation
                                temp1.alpha=log(alp)-exp(temp)*L
                                }
  
  V=temp1.alpha/rowSums(temp1.alpha)
  
  temp1.alpha=t(basis) %*% ((1-alp[,1:(kdim-1)])*(Dk[,1:(kdim-1)]+Dk[,kdim+1]*V[,1:(kdim-1)]))
  temp2.alpha=t(basis) %*% ((1-alp[,1:(kdim-1)])*alp[,1:(kdim-1)]/alp[,kdim]*(Dk[,kdim]+Dk[,kdim+1]*V[,kdim]))
  
  
  newC.coe=apply(rbind(c(temp1.alpha),c(temp2.alpha)),2,function(list) if(min(list)<0) abs(min(list))+alpha.stable.parameter+1 else alpha.stable.parameter)
  C.coe=matrix(apply(rbind(newC.coe,c(C.coe)),2,max),zdim+1,kdim-1)
  coe=coe+log(temp1.alpha+C.coe)-log(temp2.alpha+C.coe)
  
  stop1=log(temp1.alpha+C.coe)-log(temp2.alpha+C.coe)
  
######################## estimate beta ##########################
    
  temp1.beta=exp(basis %*% coe)
  temp2.beta=temp1.beta/(1+rowSums(temp1.beta))
  
  alp=cbind(temp2.beta,1-rowSums(temp2.beta))
  
  temp=Data.Z %*% b
  temp1.beta=alp*exp(-exp(temp)*L)
  
  if(sum(rowSums(temp1.beta)==0)!=0){
                                message("The value of W is zero, numerical approximation is applied.")  ## using log transformation
                                temp1.beta=log(alp)-exp(temp)*L
                                }
  
  V=temp1.beta/rowSums(temp1.beta)
  
  temp1.beta=t(Data.Z) %*% Dk[,1:kdim]
  temp2.beta=t(Data.Z) %*% (L*exp(temp)*(Dk[,1:kdim]+Dk[,kdim+1]*V[,1:kdim]))
  
  
  newC.b=apply(rbind(c(temp1.beta),c(temp2.beta)),2,function(list) if(min(list)<0) abs(min(list))+beta.stable.parameter+1 else beta.stable.parameter)
  C.b=matrix(apply(rbind(newC.b,c(C.b)),2,max),zdim,kdim)
  b=b+log(temp1.beta+C.b)-log(temp2.beta+C.b)
  
  stop2=log(temp1.beta+C.b)-log(temp2.beta+C.b)
  
######################## stopping threshold #########################

norm=c(abs(stop2),abs(c(stop1)))
if(max(norm)<=threshold){
                         message(paste("NPMLE convergence at iteration = ",i,sep=""))
                         break()
                         }
if(i==iteration)
message("attained iteration maximum")

}

Sigma.entryij=function(i,j) pLikelihood(c(b,coe)+diag(length(c(b,coe)))[i,]/sqrt(N)+diag(length(c(b,coe)))[,j]/sqrt(N),L)
Sigma.entryi=function(i) pLikelihood(c(b,coe)+diag(length(c(b,coe)))[i,]/sqrt(N),L)

V.Sigma.entryij=Vectorize(Sigma.entryij)
sigmaij=outer(seq(c(b,coe)),seq(c(b,coe)),V.Sigma.entryij)
sigmai=sapply(seq(c(b,coe)),Sigma.entryi)

sigma=-sigmaij+matrix(sigmai,length(c(b,coe)),length(c(b,coe)))+t(matrix(sigmai,length(c(b,coe)),length(c(b,coe))))-Likelihood(c(b,coe),L)

if(det(sigma)!=0 & !is.na(det(sigma))){
                  SD<-sqrt(sapply(diag(solve(sigma)),function(element) max(element,0)))/sqrt(N)
                  # SD<-sqrt(sapply(diag(solve(sigma)),function(element) abs(element)))/sqrt(N)
                  se.alpha<-matrix(SD[(kdim*zdim+1):(zdim*kdim+(zdim+1)*(kdim-1))],zdim+1,kdim-1)
                  se.beta<-matrix(SD[1:(kdim*zdim)],zdim,kdim)
                  pvalue.alpha<-1-pchisq((coe/se.alpha)^2,1)
                  pvalue.beta<-1-pchisq((b/se.beta)^2,1)}
else{
     se.alpha<- matrix(NA,zdim+1,kdim-1)
     se.beta<- matrix(NA,zdim,kdim)}

rownames(coe)=rownames(se.alpha)=rownames(pvalue.alpha)=c("(Intercept)",colnames(DATA)[3:(2+zdim)])
colnames(coe)=colnames(se.alpha)=colnames(pvalue.alpha)=paste("risk.factor",levels(factor(DATA[,2])),sep="")[1:(kdim-1)]

rownames(b)=rownames(se.beta)=rownames(pvalue.beta)=colnames(DATA)[3:(2+zdim)]
colnames(b)=colnames(se.beta)=colnames(pvalue.beta)=paste("risk.factor",levels(factor(DATA[,2])),sep="")[1:kdim]

rownames(L)<-paste("time",seq(L[,1]),sep="")
colnames(L)<-paste("risk.factor",levels(factor(DATA[,2])),sep="")[1:kdim]

list(alpha=coe,alpha.se=se.alpha,alpha.pvalue=pvalue.alpha,
     alpha.95.lower.CI=coe-qnorm(0.975)*se.alpha,alpha.95.upper.CI=coe+qnorm(0.975)*se.alpha,
     beta=b,beta.se=se.beta,beta.pvalue=pvalue.beta,
     beta.95.lower.CI=b-qnorm(0.975)*se.beta,beta.95.upper.CI=b+qnorm(0.975)*se.beta,
     Lambda=L)
}



