library(stats)


getDemand <- function(N,p,b,i){
  if(b==0)
    print("Error! b cant't be 0.")
  return (N*(1-p/b)+RandUni[i])
}


getRegret <- function(T,N,P,b,PP){
  
  #PP: Price Pulled
  
  RegVec=rep(0,T)
  for(i in 1:T){
    if(i!=1)
      RegVec[i]=RegVec[i-1]+ max(P*N*(1-P/b[i]))-(PP[i]*N*(1-PP[i]/b[i]))
    else
      RegVec[i]=max(P*N*(1-P/b[i]))-(PP[i]*N*(1-PP[i]/b[i]))
  }
  return (RegVec)
}



getQuality <- function(cq,DemObs,PP,NF=0.2,w=0.5){
  
  #cq: Cumulative quality
  #DemObs: Demand Observed in the window
  #PP: Vector of price pulled
  #NF: Normalising factor (between 0 and 1)
  mid=-1
  for(i in (1:(length(PP)-1))){
    if(PP[i]==PP[length(PP)])
    {
      mid=i
      break
    }
  }
  if(mid==-1)
  {
    return (cq)
  }
  d1=min(DemObs[length(DemObs)],DemObs[mid])
  d2=max(DemObs[length(DemObs)],DemObs[mid])
  iq = 2*exp(-(((d2-d1)/(NF*d1))^2))-1
  cq =w*cq+(1-w)*iq
  #print(iq)  
  return (cq)
}

getMeanQuality <- function(cq,DemObs,PP,NF=0.2,w=0.5,tempco,lim=10){
  
  #cq: Cumulative quality
  #DemObs: Demand Observed in the window
  #PP: Vector of price pulled
  #NF: Normalising factor (between 0 and 1)
  
  d1 = 0
  co1 = 0
  for(i in (ceiling(length(PP)/2)):1){
    if(PP[i]==PP[length(PP)])
    {
      d1=d1+DemObs[i]
      co1=co1+1
    }
  }
  
  d2 = 0
  co2 = 0
  for(i in (ceiling(length(PP)/2)):length(PP)){
    if(PP[i]==PP[length(PP)])
    {
      d2=d2+DemObs[i]
      co2=co2+1
    }
  }
  if(co1==0 || co2==0)
    return (cq)
  dd1=d1/co1
  dd2=d2/co2
  d1=min(dd1,dd2)
  d2=max(dd1,dd2)
  
  iq = (2*exp(-(((d2-d1)/(0.1*NF*d1))^2))-1)*min(1,min(co1,co2)/lim)
  cq =w*cq+(1-w)*iq
  print(c(tempco,(d2-d1),co1,co2,iq))
  return (cq)
}



sqrtMSE <- function (X,Y,fit)
{
  V2 = X
  return (sqrt(sum((predict(fit,newdata=as.data.frame(V2))-Y)^2)/(length(X)-2)))
}

findConfInt <- function (X,Y,fit,curX,alpha)
{
  N = length(X)
  if(N <= 2)
    return (Inf)
  s = sqrtMSE(X,Y,fit)
  meanX = mean(X)
  sy= s*sqrt((1/N)+(((curX-meanX)^2)/(sum((X-meanX)^2))))
  #print (c(curX,sy,s,sqrt((1/N)+(((curX-meanX)^2)/(sum((X-meanX)^2)))),((curX-meanX)^2)/(sum((X-meanX)^2)),(curX-meanX)^2,curX,meanX))
  return (qt(1-alpha/2,df=N-2)*sy)
}


T = 130 #Horizon
N = 800 #Market Size
StdDev = 40 #Standard Deviation
Tau = 20 #Sliding window length
P = c(2.0,3.0,4.0)
b = c(rep(5.5,40),rep(4.5,50),rep(9.0,40))
#b=c(rep(5.5,36),seq(5.4,4.6,-0.1),rep(4.5,40),seq(5,9,0.5),rep(9.0,36))


set.seed(1)
#RandUni = runif(T,-StdDev,StdDev) #Random error vector
RandUni = rnorm(T,0,StdDev) #Random error vector

ArmPulled = rep(0,T) #Arm pulled
PPulled = rep(0,T) #P[ArmPulled]
DemandObseved = rep(0,T) #Observed Demand
CumRew = rep(0,T) #Cumulative Reward

IQ = rep(0,T) #Vector of instantaneous quality
CQ = rep(0,T) #Vector of cumulative quality
#CQmean = rep(0,T)
##Greedy
for (i in 1:T){
  curP=0
  curT=min(Tau,i-1)
 if(i <=2){
   curP=i
 } else if (sd(ArmPulled[(i-curT):(i-1)])==0) {
    if(ArmPulled[i-1]==1){
    curP=2} else {curP=1}
  } else{
    dfr = as.data.frame(cbind(DemandObseved[(i-curT):(i-1)] , (PPulled[(i-curT):(i-1)])))
    fit = lm (V1 ~ V2 , dfr)
    V2=P
    predDem = predict(fit,as.data.frame(V2))
    #print(predDem)    
    predRev = predDem*P
    curP = which.max(predRev)
  }
  ArmPulled[i] = curP
  PPulled[i] = P[curP]
  DemandObseved[i] = getDemand(N,P[curP],b[i],i)
  CumRew[i] = CumRew[max(i-1,1)] + PPulled[i]*DemandObseved[i]
  if(i>1){
    #dfr = as.data.frame(cbind(DemandObseved[(i-curT):(i-1)] , (PPulled[(i-curT):(i-1)])))
    #NormFact=2*sd(DemandObseved[(i-curT):(i)])/mean(DemandObseved[(i-curT):(i)])
    #if(length(which(ArmPulled[i]==ArmPulled[(i-curT):(i)]))>1){
      
    #NormFact=2*sd(DemandObseved[(i-curT):(i)][ArmPulled[i]==ArmPulled[(i-curT):(i)]])/mean(DemandObseved[(i-curT):(i)][ArmPulled[i]==ArmPulled[(i-curT):(i)]])  
    NormFact=0.25  
    #print(c(NormFact,mean(DemandObseved[(i-curT):(i)][ArmPulled[i]==ArmPulled[(i-curT):(i)]])))
    IQ[i]=getQuality(CQ[i-1],DemandObseved[(i-curT):(i)],PPulled[(i-curT):(i)],NormFact,w=0.5)
    #CQmean[i]=getMeanQuality(CQmean[i-1],DemandObseved[(i-curT):(i)],PPulled[(i-curT):(i)],NormFact,w=0.5)
    #}
    
  }
  else{
    IQ[i]=1
    #CQmean[i]=1
  }
  
}
RegVec = getRegret(T,N,P,b,PPulled)
plot(c(1:130),RegVec,type="l",col="red",xlab="t",ylab="Regret")
##Greedy

##Greedy with window cutting
Tau=0
ArmPulled = rep(0,T) #Arm pulled
PPulled = rep(0,T) #P[ArmPulled]
DemandObseved = rep(0,T) #Observed Demand
CumRew = rep(0,T) #Cumulative Reward
TauRec = rep(0,T)
IQ = rep(0,T) #Vector of instantaneous quality
CQ = rep(0,T) #Vector of cumulative quality
for (i in 1:T){
  curP=0
  curT=Tau
  TauRec[i]=Tau
  if(Tau < 2){
    if(Tau==0){
    curP=Tau+1}
    else{
      curP=((ArmPulled[i-1]+1)%%(length(P)))+1}
  } else if (sd(ArmPulled[(i-curT):(i-1)])==0) {
    if(ArmPulled[i-1]==1){
      curP=2} else {curP=1}
  } else{
    dfr = as.data.frame(cbind(DemandObseved[(i-curT):(i-1)] , (PPulled[(i-curT):(i-1)])))
    fit = lm (V1 ~ V2 , data=dfr)
    
    V2=P
    predDem = predict(fit,newdata=as.data.frame(V2))
    #print(predDem)    
    predRev = predDem*P
    curP = which.max(predRev)
    for (coi in 1:length(P))
    {
    print(c(Tau, P[coi]*findConfInt((PPulled[(i-curT):(i-1)]),DemandObseved[(i-curT):(i-1)] ,fit,P[coi],0.05), predRev[coi], sqrtMSE(PPulled[(i-curT):(i-1)],DemandObseved[(i-curT):(i-1)] ,fit)))
    }
  }
  ArmPulled[i] = curP
  PPulled[i] = P[curP]
  DemandObseved[i] = getDemand(N,P[curP],b[i],i)
  CumRew[i] = CumRew[max(i-1,1)] + PPulled[i]*DemandObseved[i]
  if(i>1){
    #NormFact=40/400
    NormFact=2*sd(DemandObseved[(i-curT):(i)])/mean(DemandObseved[(i-curT):(i)])
    CQ[i]=getQuality(CQ[i-1],DemandObseved[(i-curT):(i)],PPulled[(i-curT):(i)],NormFact,w=0.5)
    if(CQ[i]<0)
    {
      Tau=ceiling(Tau/2)
      #CQ[i]=0
    }
    #print (Tau)
  }
  else{
    CQ[i]=1
    #CQmean[i]=1
  }
  Tau=Tau+1
}

RegVec = getRegret(T,N,P,b,PPulled)
lines(c(1:130),RegVec,type="l",col="green")
