library(MASS)
library(splines)

set.seed(10)
p<-2
N<-100

sigma<-diag(rep(1,p))
mu<-rep(0,p)

X<-mvrnorm(N,mu,sigma)
Y<-rnorm(N,mean=0,sd=5)


ppr_function<-function(X,Y){
  d<-1
  p<-ncol(X)
  N<-nrow(X)
  w<-rep(1/p,p)
  w<-w/sqrt(sum(w^2))
  
  while(d>=10^(-5)){
    v<-X%*%w
    v_order<-sort(v)
    cubic_base<-matrix(nrow=N,ncol=N+4)
    cubic_base[,1]<-1
    cubic_base[,2]<-v
    cubic_base[,3]<-v^2
    cubic_base[,4]<-v^3
    for(i in 1:N){
      cubic_base[,4+i]<-apply(v,1,function(x)ifelse((x>v[i]),1,0))*(v-v_order[i])^3
    }
    
    X_N<-matrix(nrow=N,ncol=N)
    X_N[,1]<-1
    X_N[,2]<-v
    for(i in 1:(N-2)){
      X_N[,i+2]<-(cubic_base[,i+4]-cubic_base[,N+4])/(v_order[i]-v_order[N])-(cubic_base[,N+3]-cubic_base[,N+4])/(v_order[N-1]-v_order[N])
    }
    
    omega<- function(x) {
      N<-length(x)
      Omega <- matrix(rep(0, N*N), ncol = N)
      knot<-sort(x)
      Omega.tmp<-matrix(rep(0, (N-2)*(N-2)), ncol = N-2)
      
      for (k in 1:(N-2)) {
        for (j in k:(N-2)) {
          K<-length(knot)
         
          C_k1 <- 1/(knot[k]-knot[K])
          C_k2 <- 1/(knot[K-1]-knot[K])
          C_j1 <- 1/(knot[j]-knot[K])
          C_j2 <- 1/(knot[K-1]-knot[K])
          integrat1<-function(x) {36*C_k1*(x-knot[k])*C_j1*(x-knot[j])}
          v1 <- integrate(integrat1, lower = knot[j], upper=knot[K-1])$value
          integrat2<- function(x) {(6*C_k1*(x-knot[k])- 6*C_k2*(x-knot[K-1]))*(6*C_j1*(x-knot[j])- 6*C_j2*(x-knot[K-1]))}
          v2<-integrate(integrat2, lower = knot[K-1], upper=knot[K])$value
          Omega.tmp[k, j] <- v1+v2
          Omega.tmp[j, k] <- Omega.tmp[k, j]
        }
      }
      Omega[3:N,3:N]<-Omega.tmp
      return(Omega)
    }
    
    lambda<-10
    Omega<-omega(v)
    beta<-solve(t(X_N)%*%X_N+lambda*Omega)%*%t(X_N) %*%Y
    
    y_hat<-X_N%*%beta
    
    
    ###calculate g'(v)
    de_mat<-matrix(nrow=N,ncol=N-2)

    for(i in 1:N){
      for(j in 1:(N-2)){
        if(j<=rank(v)[i]){
          de_mat[i,j]<-3*(v[i]-v_s[j])^2/(v_s[j]-v_s[N])
        }else{
          de_mat[i,j]<-0
        }
      }
    }

    de<-matrix(nrow=N,ncol=1)
    de<-matrix(beta[2],nrow=N,ncol=1)+de_mat%*%beta[3:N,1]
    
    s<-v+(Y-y_hat)/de
    
    W<-diag(c(de))
    
    w_new<-ginv(t(X)%*%W%*%X)%*%t(X)%*%W%*%s

    d<-sqrt(sum((w-w_new)^2))
    w<-w_new
    w<-w/sqrt(sum(w^2))
  }
  print(w)
}
ppr_function(X,Y)

f1<-ppr(X,Y,nterms=2)
f1$alpha[,1]
