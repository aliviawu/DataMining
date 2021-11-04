
library(mvtnorm)

######### Perceptron model with Gradient Descent ###########
#GD
w_old<- c(1,1,1)
rho=0.01
counter<- 0 
seperate_criterion<- -y*(obs%*%w_old)
misc_set<- which(seperate_criterion>0)
w_new<- NULL
converge<- NULL

b<- 1
for(i in 1:10000){
  n<- 10000
  p<- 2
  x<- rmvnorm(n,mean=rep(0,p),sigma=(diag(1,p,p)+(1-diag(1,p,p))*0.15))
  z<- rmvnorm(n, mean=rep(6, p), sigma = (diag(1,p,p)+(1-diag(1,p,p))*0.1)) 
  obs<- rbind(x, z)
  y<- rep(c(-1,1), each=n)
  obs<- cbind(1,obs)
  w_old<- c(1,1,1)
  rho=0.01
  seperate_criterion<- -y*(obs%*%w_old)
  misc_set<- which(seperate_criterion>0)
  w_new<- NULL
  b<- 1
  repeat{
    
    delta<- t(y[misc_set])%*%obs[misc_set,]
    w_new<- w_old +rho*delta
    
    seperate_criterion<- -y*(obs%*%t(w_new))
    misc_set<- which(seperate_criterion>0)
    if(length(misc_set)==0){
      converge[i]=1
      break}
    w_old=w_new
    b<- b+1
    if(b>10000){break}}}
 
sum(!is.na(converge))/10000 # 0.8602

#I tried 20000 sample size and 10000 times of simulation.The converge rate is 0.8602.


############ Perceptron model with stochastic gradient descent ##############

w_old<- c(1,1,1)
rho=0.01
counter<- 0 
seperate_criterion<- -y*(obs%*%w_old)
misc_set<- which(seperate_criterion>0)
w_new<- NULL
converge<- NULL

b<- 1
for(i in 1:10000){
  n<- 10000
  p<- 2
  x<- rmvnorm(n,mean=rep(0,p),sigma=(diag(1,p,p)+(1-diag(1,p,p))*0.15))
  z<- rmvnorm(n, mean=rep(6, p), sigma = (diag(1,p,p)+(1-diag(1,p,p))*0.1)) 
  obs<- rbind(x, z)
  y<- rep(c(-1,1), each=n)
  obs<- cbind(1,obs)
  w_old<- c(1,1,1)
  rho=0.01
  seperate_criterion<- -y*(obs%*%w_old)
  misc_set<- which(seperate_criterion>0)
  w_new<- NULL
  b<- 1
  repeat{
    
    index<- sample(misc_set, 1)
    delta<- y[index]*obs[index,]
    w_new<- w_old +rho*delta
    
    seperate_criterion<- -y*(obs%*%(w_new))
    misc_set<- which(seperate_criterion>0)
    if(length(misc_set)==0){
      converge[i]=1
      break}
    w_old=w_new
    b<- b+1
    if(b>10000){break}
   
  }}
#958/1000
#4760/5000
sum(!is.na(converge))/10000 #0.8511


#Sample size: 20000.number of simulations: 10000. The accuracy rate is 0.8511.
#Concluded, the converge rate of Gradient descent is slightly greater than stochastic gradient descent. 



