set.seed(100)
N   <- 100
Z1  <- rnorm(N, 0,1)
px1 <- 0.3
X1  <- rbinom(N, 1, px1) ; mean(X1)

beta0  <- c(-3,1,1)  #beta <- 1
#X2    <- cbind(Z1,X1)
XX    <- cbind(1, Z1, X1) # used for X in the function, gives N*3 matrix
eta   <- XX%*%beta0
mu    <- function(eta) as.vector(exp(eta)/(1+exp(eta))) # 1/(1+exp(-eta)) 

logit <-  function(p) log(p/(1-p))
 


delta <- c(1,-1,-1) #create a empty matrix for delta
Y     <-  rbinom(N, 1, mu(eta))  

## for the pi probability, if we need to assign more
###samples, then we need to change the delta_0
pi    <-  1/(1+exp(-1+Z1+ Y))  # delta= c(1,-1,-1)
sum(pi)

###take the sample from the data
Ri     <- as.numeric(runif(N, 0,1)  <= pi )
sum(Ri[Ri==1])  ### size of sample size, (small n)

element <- diag(Ri)[1:10,1:10] 

#######
DataF  <- data.frame(X1, Z1, Y, Ri, pi, wi= 1/pi)
library(dplyr)
DataS  <-  filter(DataF, Ri==1) #sampled data 

### Fit model
beta<- NR()






##############each small function for dervi###

calc_pi = function(Y,Z,delta)
{
  int <-  cbind(1,Y,Z)
  pivalue <-  as.vector(exp(int%*%delta) / (1+ exp(int%*%delta)))
  return(pivalue)
} 


derv_pi1 =function(Y,Z,delta)
{
  
  ZZ <- cbind(1,Y,Z)
  exp.part <- as.vector(exp(-int%*%delta))
  value_pi <- (exp.part * int)/ ((1+exp.part)^2)
  return(value_pi)
}

####return the value for (intercept, Y,Z) matrix Nx3 ? May I use as.vector for the poutput?
derv_delta=function(Y,Z,delta){
  ZZ <- cbind(1,Y,Z)
  
  zz.list <- split(ZZ,seq(nrow(ZZ)))
  
  Dvalue <- lapply(1:nrow(ZZ), function(i) ## Then I apply each row of (1,Yi,Zi) 
    ## in to calc_pi(Y,Z,delta) function 
  {
    
    pi1 <- calc_pi(zz.list[[i]][2],zz.list[[i]][3],delta)
    const <- pi1*(1-pi1)
    value <- const*zz.list[[i]]
    return(value)
  })
  #Dvalue <- calc_pi(Y,Z,delta)*(1-calc_pi(Y,Z,delta))*ZZ
  
  return(Dvalue)
  
}




derv_delta2 <- function(Y,Z,delta){
  #Z=Z1
  #
  
  ZZ <- cbind(1,Y,Z) ### create a big matrix (Nx3) 
  
  zz.list <- split(ZZ,seq(nrow(ZZ))) ##for each row of matrix, I make it as a list
  
  value <- lapply(1:nrow(int), function(i) ## Then I apply each row of (1,Yi,Zi) 
    ## in to calc_pi(Y,Z,delta) function 
  {
    
    pi1 <- calc_pi(zz.list[[i]][2],zz.list[[i]][3],delta)
    const <- pi1*(1-pi1)*(1-2*pi1)
    value <- const*zz.list[[i]]%*% t(zz.list[[i]]) 
    return(value)
  })
  
  #value1 <- matrix(unlist(value),ncol=3,byrow=TRUE)
  #value1 <- Reduce('+',value) ## add all element of list together
  return(value)  ###The value which return the 300x3 matrix 
  
  
}



####der_beta output: vector
der_beta <- function(X,Z,beta){
  #Z=Z1
  #X=X1
  
  int <- cbind(1,Z,X)
  
  N <- nrow(int)
  
  int.list <- split(int,seq(nrow(int)))
  
  
  Bvalue <- lapply(1:nrow(int), function(i)
  {
    
    eta <- int.list[[i]]%*%beta
    
    const <-  mu(eta)*(1-mu(eta))
    
    value <- const*int.list[[i]]
    
    return(value)
  })
  
   #value <- mu(eta)*(1-mu(eta))*int   ###matrix Nx3
 return(Bvalue) 
}

der_beta2 <- function(X,Z,beta){
  #Z=Z1
  #X=X1
  #beta=beta
  
  int <- cbind(1,X,Z)
  
  int.list <- split(int,seq(nrow(int)))
  
  value2 <- lapply(1:nrow(int), function(i)
  {
    
    eta <- int.list[[i]]%*%beta
    
    const <-  mu(eta)*(1-mu(eta))*(1-2*mu(eta))
    
    value2 <- const*int.list[[i]]%*% t(int.list[[i]]) 
    
    return(value2)
  })
  
  return(value2)
}



U  <- function(Y, X,Z, R, beta, delta){
  #Z=Z1
  #X=X1
  #R=Ri
  int <- cbind(1,X,Z)
  
  int.list <- split(int,seq(nrow(int)))
  
  eta1 <- as.vector(int%*%beta)
  piy0 <- calc_pi(0,Z,delta)
  piy1 <- calc_pi(1,Z,delta)
  
  numerator1 <- mu(eta1) *(piy0 - Y*piy1+ piy1) - piy0*Y 
  
  denominator1 <- mu(eta1)*(1- mu(eta1)) *(piy0 -  mu(eta1)*piy1 - piy0)  
  u11  <- R *(numerator1/denominator1)
  aux  <- -der_beta(X,Z,beta)
  
  Uvalue1 <- Map ('*', u11, aux)
  
  frac <- lapply(1:nrow(int), function(i){
    
  eta <- int.list[[i]]%*%beta
    
   numerator <- mu(eta)* (calc_pi(0,int.list[[i]][3],delta)- Y[i]*calc_pi(1,int.list[[i]][3],delta)
                  +calc_pi(1,int.list[[i]][3],delta))-calc_pi(0,int.list[[i]][3],delta)*Y[i]
   
   denominator <- mu(eta)*(1-mu(eta))*((calc_pi(0,int.list[[i]][3],delta)-calc_pi(1,int.list[[i]][3],delta))
                  *mu(eta) -calc_pi(0,int.list[[i]][3],delta)) 
   
   u1 <- R[i] *(numerator/denominator)
   
   Uvalue <- u1[i]*der_beta(X,Z,beta)[[i]]
   
   #return(Uvalue)
  })
  
  fracU <- Reduce('+',frac)
  ## 
  ##if(link=='exponential'){
  ##if(link=='probit'){}
 
  return(fracU)
}


####whether the output of these functions are a column of number
Uphi  <- function(Y, Z, R,  delta  ){
  #Z=Z1
  #R=Ri
  ratio <- (R-calc_pi(Y,Z,delta))/
       (calc_pi(Y,Z,delta)*(1-calc_pi(Y,Z,delta)))  ## a vector of constant
  
  der.value <- Map("*", ratio, derv_delta(Y,Z,delta))
  
  phi <-  Reduce('+',der.value)
  
  return(phi)
}



S  <- function(Y, X,Z, R, beta, delta ){
 #X <- X1
 #Z=Z1
 #R= Ri
 U1<-  U(Y, X,Z, R, beta, delta) 
 
 U2 <-   Uphi(Y, Z, R,  delta  )
 
 Svalu<-  c(U1 , U2)
 
 return(Svalu)
 
}

 
################################## Score function 
dS_beta  <- function( Y, X,Z, R, beta, delta ){
  Z=Z1
  beta=beta0
  
  int = cbind(1,X,Z)
 
  D <- matrix(unlist(der_beta(X,Z,beta)),ncol=3,byrow=T)
  
  eta <- int%*%beta
  
  expect_pi <- calc_pi(1,Z,delta)*mu(eta)   ###a 1:100 vector
             + calc_pi(0,Z,delta)*(1-mu(eta))
  
  diag1 <- as.vector(R/(mu(eta)*(1-mu(eta))))
  
  diag2 <- as.vector((R*(Y-mu(eta)))/((mu(eta)^2)*(1-mu(eta))^2))
  
  diag3 <- as.vector(1-2*mu(eta))
  
  
  exp1 <- as.vector(mu(eta)*calc_pi(1,Z,delta))*D
         +as.vector((1-mu(eta))*calc_pi(0,Z,delta))*D
  
  numerator <- mu(eta)* (calc_pi(0,Z,delta)- Y*calc_pi(1,Z,delta)+calc_pi(1,Z,delta))
             -calc_pi(0,Z,delta)*Y  ###a 1:100 vector
  
  denominator <- mu(eta)*(1-mu(eta))*((calc_pi(0,Z,delta)-calc_pi(1,Z,delta))*mu(eta) 
                -calc_pi(0,Z,delta)) ###a 1:100 vector
   
  aux0 <- R *(numerator/denominator)
  aux <- der_beta2(X,Z,beta)
                               
  u1 <- Map('*', aux0, aux )
  sum_u1 <- Reduce('+',u1)  ####colSums(u1)
  
  Svalue <- -t(D*(diag1+diag2*diag3))%*%D+
            t(D)%*%(R*(1/expect_pi)*exp1) +sum_u1  ##matrix 3x3
  
  return(Svalue)
} 




dS_delta <- function(X,Z,R,beta,delta){
  #Z=Z1
  #beta=beta0
  
  int <- cbind(1,X,Z)
  
  D <- matrix(unlist(der_beta(X,Z,beta)),ncol=3,byrow=T)
  
  eta <- int%*%beta
  
  expD_tilde <- lapply(1:N, function(i)
    { derv_delta(1,Z,delta)[[i]]+derv_delta(0,Z,delta)[[i]] 
      })##matrix Nx3
  
 D_tilde <- matrix(unlist(expD_tilde),ncol = 3,byrow =TRUE )
  
  
  expect_pi <- calc_pi(1,Z,delta)*mu(eta) 
                + calc_pi(0,Z,delta)*(1-mu(eta))  ##A vector
  
  frac1 <- diag((1-mu(eta))/(mu(eta)*(1-mu(eta))))  ##matrix NxN
  
  frac0 <- diag((0-mu(eta))/(mu(eta)*(0-mu(eta))))   ##matrix NxN
  
  #ddpart <- matrix(unlist(derv_delta()))
  
  exp1 <- frac1%*%matrix(unlist(derv_delta(1,Z,delta)),ncol = 3,byrow = TRUE)*mu(eta)
        + frac0%*%matrix(unlist(derv_delta(0,Z,delta)),ncol = 3,byrow = TRUE)*(1-mu(eta))   ##matrix Nx3
  
  exp2 <- frac1%*%diag(calc_pi(1,Z,delta)*mu(eta))
        + frac0%*%diag(calc_pi(0,Z,delta)*(1-mu(eta))) ##matrix 100x100
  
  value <- -t(D)%*%diag(1/expect_pi)%*%diag(R)%*%exp1
          + t(D)%*%diag(1/expect_pi^2)%*%diag(R)%*%exp2%*%D_tilde ##matrix 3x3
  
  return(value)
}


dphi_delta <- function(R,Y,Z,delta){
  #R=Ri
  #Z=Z1
  
  numerator <- R-calc_pi(Y,Z,delta) ## a vector 1:100
  
  denominator <- calc_pi(Y,Z,delta)*(1-calc_pi(Y,Z,delta)) ## a vector 1:100
  
  part1 <- numerator/denominator #a vector 1:100
  
  part2 <- (-1/denominator)- ((numerator*(1-2*calc_pi(Y,Z,delta)))/denominator^2)#a vector 1:100
  
  v.derpi2 <- lapply(1:N, function(i){
    
    value.part1 <- part1[i]*derv_delta2(Y,Z,delta)[[i]]
    
    return(value.part1 )
  }) # a 3x3 matrix for N
  
  
  v.derpi1 <- lapply(1:N, function(i){
    
    value.part2 <- part2[i]*derv_delta(Y,Z,delta)[[i]]%*%t(derv_delta(Y,Z,delta)[[i]])
    
    return(value.part2 )
    })## matrix 100x100 wrong fix
    
   

  Phivalue <- Reduce("+",v.derpi1) + Reduce("+",v.derpi2) ## matrix100x100
  
  return(Phivalue)
}





################Newton-Raphson method########

#Inputs:
set.seed(1)
N   <- 1000
Z1  <- rnorm(N, 0,1)
px1 <- 0.3
X1  <- rbinom(N, 1, px1)
XX    <- cbind(1, Z1, X1)
beta0  <- c(1,1,1 )
eta0   <- XX%*%beta0
mu0    <- function(eta0) as.vector(exp(eta0)/(1+exp(eta0))) # 1/(1+exp(-eta)) 

Ri <- as.numeric(runif(N, 0,1)  <= pi )
delta <- c(1,-1,-1) #create a empty matrix for delta
Y     <-  rbinom(N, 1, mu(eta0))  

#Initial values
delta <- rep(0.1, 3)
beta <- glm(Y~X1+Z1, family = binomial)$coefficients
epsi <- 1
eps0  <- 0.000001

theta_n <-  c(beta, delta)

counti  <- 1
#Newton-Raphson method:

while(epsi  > eps0  )
{
  f  <-  S(Y, X,Z, R, beta, delta )   ## S_C, psi
  
  
  #Derivative of dS_epslon w.r.t. beta:
  d11 <- dS_beta(Y,X,Z,R,beta,delta)  ## d S_C/ d beta
  d12 <- dS_delta(X,Z,R,beta,delta)  ## d S_C/ d delta
  
  
  d21  <- matrix(0,ncol = 3,nrow = 3,byrow = TRUE) ## d psi/ d beta, whether it is a matrix
  d22  <- dphi_delta(R,Y,Z,delta)  ## d psi/ d delta
    
    #Derivative of dS_epsilon w.r.t. delta: 
    
    
    #Derivative of f(eps):
    df  <-   rbind( cbind(d11, d12), cbind(d21, d22)    )
  
  theta_n1  <- theta_n  -     solve(df) %*%  f
  
  
  #Update epsilon: 
  epsi  <- max(abs(theta_n1 - theta_n ))
  
  beta  <-  theta_n1[1:3]
  delta  <-  theta_n1[4:6]
  
  print(counti)
  print(epsi)
  print(theta_n1)
  counti  <- counti +1
  
}

