
set.seed(1)
T = 500 #Horizon
N = 800 #Market Size
#StdDev randomly between 5 and 50
#gamma between 1 and 20
#gamma points between 2 and 499
#gamma + 1 b values

for (i in 1:20){
  
  StdDev = runif(1,5,51)
  gamma = ceiling(runif(1,0,20))
  cpoints = sort(sample(2:500)[1:gamma])
  bvals = runif(gamma+1,4.5,10)
  print(c("T=",T))
  print(c("N=",N))
  print(c("StdDev=",StdDev))
  print(c("NumBreaks=",gamma))
  print("c(")
  for (j in 1:(gamma+1)){
    if(j==1){
      print(c("rep(",bvals[j],",",cpoints[j]-1,"),"))  
    }
    else if(j==gamma+1)
    {
      print(c("rep(",bvals[j],",",T-cpoints[j-1]+1,"),"))  
    }
    else
    {
    print(c("rep(",bvals[j],",",cpoints[j]-cpoints[j-1],"),"))
    }
  }
  print(")")
  print(c("BreakPoints=",cpoints))
  print(c("bvals = ",bvals))
  print("-------------------------------------------")
  
}