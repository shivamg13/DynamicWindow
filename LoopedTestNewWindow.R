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
  #print(c(tempco,(d2-d1),co1,co2,iq))
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
  
  s = sqrtMSE(X,Y,fit)
  meanX = mean(X)
  sy= s*sqrt((1/N)+(((curX-meanX)^2)/(sum((X-meanX)^2))))
  #print (c(curX,sy,s,sqrt((1/N)+(((curX-meanX)^2)/(sum((X-meanX)^2)))),((curX-meanX)^2)/(sum((X-meanX)^2)),(curX-meanX)^2,curX,meanX))
  if(N <= 2 || is.infinite(sy) || is.nan(sy))
  {
    return (Inf)
  }
  return (qt(1-alpha/2,df=N-2)*sy)
}


findConfIntPoint <- function (X,Y,fit,curX,alpha)
{
  N = length(X)
  
  s = sqrtMSE(X,Y,fit)
  meanX = mean(X)
  sy= s*sqrt(1+(1/N)+(((curX-meanX)^2)/(sum((X-meanX)^2))))
  #print (c(curX,sy,s,sqrt((1/N)+(((curX-meanX)^2)/(sum((X-meanX)^2)))),((curX-meanX)^2)/(sum((X-meanX)T^2)),(curX-meanX)^2,curX,meanX))
  if(N <= 2 || is.infinite(sy) || is.nan(sy))
  {return (Inf)}
  return (qt(1-alpha/2,df=N-2)*sy)
}


winInc <- function(X,Y ,fit,Dom,alpha,tHold){
  Z = findConfInt(X,Y,fit,Dom,alpha)
  
  if(max(Z)==Inf)
  {
    return (TRUE)
  }
  Z = Z*Dom
  V2 = Dom
  Pred =  Dom*predict(fit,newdata=as.data.frame(V2))
  U = (Pred + Z - max(Pred))/max(Pred)
  #print(Pred)
  #print(Z)
  if(max(U) > tHold)
  {
    return (TRUE)
  }
  else
  {
    return (FALSE) 
  }
}


findSdConfInt <- function (X,Y,fit,alpha)
{
  V2 = X
  lenX = length(X)
  s2 = (sum((predict(fit,newdata=as.data.frame(V2))-Y)^2))
  u_int = sqrt(s2/(qchisq(alpha/2, lenX -2 )))
  l_int = sqrt(s2/(qchisq(1-alpha/2, lenX -2 )))
  return (c(l_int,u_int))
}


getConfQuality <- function(prevCQ, Y, X, alpha1, wt,iter,TauVal)
{
  
  mid = (floor((1+length(X))/2))
  #print(c(i, Tau , sqrtMSE((PPulled[(i-curT):(i)]), DemandObseved[(i-curT):(i)],fit),sqrtMSE((PPulled[(i-curT):(mid)]),DemandObseved[(i-curT):(mid)] ,fit),sqrtMSE((PPulled[(mid+1):(i)]),DemandObseved[(mid+1):(i)] ,fit)))
  #print (c(i,Tau,findSdConfInt((PPulled[(i-curT):(mid)]),DemandObseved[(i-curT):(mid)] ,fit,0.25),findSdConfInt((PPulled[(mid+1):(i)]),DemandObseved[(mid+1):(i)] ,fit,0.25)))
  dfr1 = as.data.frame(cbind(Y[1:mid] , X[1:mid]))
  
  
  if((min(dfr1$V2)==max(dfr1$V2)))
  {
    if((min(dfr1$V2)==X[length(X)]))
    {
      predDem1=mean(dfr1$V1)
      N1 = length(dfr1$V2)
      CInt = (sd(dfr1$V1)*qt(1-alpha1/2,df=N1-1))/sqrt(N1)
    }
    else
    {
      predDem1 = Inf
      CInt = Inf
    }
    
  }
  else
  {
    fit1 = lm (V1 ~ V2 , data=dfr1)
    V2=X[length(X)]
    predDem1 = predict(fit1,newdata=as.data.frame(V2))
    CInt = findConfIntPoint((X[(1):(mid)]),Y[(1):(mid)],fit1,V2,alpha1)
  }
  if(is.infinite(predDem1) || is.nan(predDem1) || is.infinite(CInt) || is.nan(CInt))
  {
    #print(c(iter, TauVal, predDem1,CInt,DemandObseved[i], Inf))
    return (prevCQ)
  }
  ratio1 = abs((Y[length(Y)] - predDem1)/(CInt))*0.6931472  #Multiplied by -log(0.5)
  IQ = 2*exp(-(ratio1)^2) - 1
  #print(c(iter,ratio1,IQ,mid))
  ##print(c(iter, TauVal, predDem1,CInt,DemandObseved[i], IQ))
  return (wt*prevCQ + (1-wt)*IQ)
  
  
}

getNewConfQuality <- function(prevCQ, Y, X, alpha1)
{
  
  mid = (floor((1+length(X))/2))
  #print(c(i, Tau , sqrtMSE((PPulled[(i-curT):(i)]), DemandObseved[(i-curT):(i)],fit),sqrtMSE((PPulled[(i-curT):(mid)]),DemandObseved[(i-curT):(mid)] ,fit),sqrtMSE((PPulled[(mid+1):(i)]),DemandObseved[(mid+1):(i)] ,fit)))
  #print (c(i,Tau,findSdConfInt((PPulled[(i-curT):(mid)]),DemandObseved[(i-curT):(mid)] ,fit,0.25),findSdConfInt((PPulled[(mid+1):(i)]),DemandObseved[(mid+1):(i)] ,fit,0.25)))
  dfr1 = as.data.frame(cbind(Y[1:mid] , X[1:mid]))
  
  
  if((min(dfr1$V2)==max(dfr1$V2)))
  {
    if((min(dfr1$V2)==X[length(X)]))
    {
      predDem1=mean(dfr1$V1)
      N1 = length(dfr1$V2)
      CInt = (sd(dfr1$V1)*qt(1-alpha1/2,df=N1-1))/sqrt(N1)
    }
    else
    {
      predDem1 = Inf
      CInt = Inf
    }
    
  }
  else
  {
    fit1 = lm (V1 ~ V2 , data=dfr1)
    V2=X[length(X)]
    predDem1 = predict(fit1,newdata=as.data.frame(V2))
    CInt = findConfIntPoint((X[(1):(mid)]),Y[(1):(mid)],fit1,V2,alpha1)
  }
  if(is.infinite(predDem1) || is.nan(predDem1) || is.infinite(CInt) || is.nan(CInt))
  {
    #print(c(iter, TauVal, predDem1,CInt,DemandObseved[i], Inf))
    return (1)
  }
  ratio1 = abs((Y[length(Y)] - predDem1)/(CInt))  #Multiplied by -log(0.5)
  RVar=1
  if(ratio1>1){
    RVar=0
  }
  #print(ratio1)
  return (RVar)
  #print(c(iter,ratio1,IQ,mid))
  ##print(c(iter, TauVal, predDem1,CInt,DemandObseved[i], IQ))
  #return (wt*prevCQ + (1-wt)*IQ)
  
  
}


T = 500 #Horizon
N = 800 #Market Size
P = c(2.0,3.0,4.0)
StdDev=          44.9152282793075
NumBreaks= 6         
b=c(
  rep(             8.34609719563741 ,                89               ),              
  rep(             9.44672616838943 ,                28               ),              
  rep(             7.68506058887579 ,                8                ),              
  rep(             8.78601149783935 ,                59               ),              
  rep(             6.54962118004914 ,                11               ),              
  rep(             7.97926613478921 ,                172              ),              
  rep(             5.69709108141251 ,                133              )
)
for (TauVal in c(seq(10,50,10),75,100)){
  #for (wtVal in seq(0.75,0.75,0.1)){
    CurReg = 0
    CurRegGreedy = 0
    for(seedVal in 1:10){


Tau = TauVal #Sliding window length

#b = c(rep(5.5,40),rep(4.5,50),rep(9.0,40))
#b=c(rep(5.5,36),seq(5.4,4.6,-0.1),rep(4.5,40),seq(5,9,0.5),rep(9.0,36))


set.seed(seedVal)
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
  
  
}
RegVec = getRegret(T,N,P,b,PPulled)
#plot(c(1:130),RegVec,type="l",col="red",xlab="t",ylab="Regret",ylim=c(0,14000))
CurRegGreedy = CurRegGreedy + RegVec[T]

    }
    print(c(TauVal,CurRegGreedy/10.0))
    }

##Greedy

for (confVal in c(0.01,0.05,0.1,0.2,0.25)){
  #for (wtVal in seq(0.75,0.75,0.1)){
  CurReg = 0
  CurRegGreedy = 0
  for(seedVal in 1:10){


    set.seed(seedVal)
    #RandUni = runif(T,-StdDev,StdDev) #Random error vector
    RandUni = rnorm(T,0,StdDev) #Random error vector
    
##Greedy with window cutting
Tau=0
ArmPulled = rep(0,T) #Arm pulled
PPulled = rep(0,T) #P[ArmPulled]
DemandObseved = rep(0,T) #Observed Demand
CumRew = rep(0,T) #Cumulative Reward
TauRec = rep(0,T)
RVarRec = rep(1,T)
IQ = rep(0,T) #Vector of instantaneous quality
CQ = rep(0,T) #Vector of cumulative quality
lastReset = 0
for (i in 1:T){
  curP=0
  curT=Tau
  TauRec[i]=Tau
  if(Tau < 2){
    if(Tau==0){
      curP=Tau+1}
    else{
      curP=((ArmPulled[i-1]+1)%%(length(P)))+1}
  } else if (min(ArmPulled[(i-curT):(i-1)])==max(ArmPulled[(i-curT):(i-1)])) {
    #if(abs(ArmPulled[i-1]-1) < abs(ArmPulled[i-1]-length(P))){
    #  curP=length(P)} else {curP=1}
    if(ArmPulled[i-1] == 1)
    {
      curP = 2
    }
    else
    {
      curP = ArmPulled[i-1] - 1 
    }
  } else{
    dfr = as.data.frame(cbind(DemandObseved[(i-curT):(i-1)] , (PPulled[(i-curT):(i-1)])))
    fit = lm (V1 ~ V2 , data=dfr)
    
    V2=P
    predDem = predict(fit,newdata=as.data.frame(V2))
    #print(predDem)    
    predRev = predDem*P
    curP = which.max(predRev)
    x=5
    #print(c(Tau, P*findConfInt((PPulled[(i-curT):(i-1)]),DemandObseved[(i-curT):(i-1)] ,fit,P,0.05), predRev, sqrtMSE(PPulled[(i-curT):(i-1)],DemandObseved[(i-curT):(i-1)] ,fit)))
    
  }
  
  ArmPulled[i] = curP
  PPulled[i] = P[curP]
  DemandObseved[i] = getDemand(N,P[curP],b[i],i)
  CumRew[i] = CumRew[max(i-1,1)] + PPulled[i]*DemandObseved[i]
  ##if(i>1){
  #NormFact=40/400
  if(Tau<3)
  {
    Tau = Tau + 1
    CQ[i] = 1
  }
  else{
    #NormFact=2*sd(DemandObseved[(i-curT):(i)])/mean(DemandObseved[(i-curT):(i)])
    #CQ[i]=getQuality(CQ[i-1],DemandObseved[(i-curT):(i)],PPulled[(i-curT):(i)],NormFact,w=0.5)
    CQ[i] = getConfQuality(CQ[i-1],DemandObseved[(i-curT):(i)],PPulled[(i-curT):(i)],0.1,w=0.75,i, Tau)
    
    ###Other method    
    #alphaR = 0.4  
    #RVarRec[i]=  getNewConfQuality(CQ[i-1],DemandObseved[(i-curT):(i)],PPulled[(i-curT):(i)],alphaR)
    #RVarWidth = 10
    #if(i-lastReset < RVarWidth)
    #{
    #  CQ[i]=1
    #}
    #else
    #{
    # CQ[i]=mean(RVarRec[(i-RVarWidth+1):i])
    #print(c(i,RVarRec[i],CQ[i],curT))
    #}
    
    #if(CQ[i]<1-alphaR-0.1)
    #{
    #  Tau = 1
    #  lastReset=i
    #}
    ######Other method ends
    if(CQ[i]<0)
    {
      Tau=1
      CQ[i]=1
    }
    else if (winInc((PPulled[(i-curT):(i-1)]),DemandObseved[(i-curT):(i-1)] ,fit,P,0.05,confVal) == TRUE)
    {
      Tau=Tau+1
    }
  }
  
  
  
}

RegVec = getRegret(T,N,P,b,PPulled)
CurReg = CurReg + RegVec[T]
    }
    print(c(confVal,CurReg/10.0))
    #print(CurReg/10.0)
  }
  
#lines(c(1:130),RegVec,type="l",col="green")
