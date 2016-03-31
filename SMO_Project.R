# Author:       Uwe Hoenig / Team DUN-Y
# Course:       SMO
# Last update:  31.03.16
# Type:         Project

#In this document I have listed a few functions that compute the price of an 
#(US/EU) Call Option

###PARAMETERS
#For simplicity let us assume no Dividends are being paid, even though it can
#be adjusted for that

#R Codes for Black-Scholes Formula to Price European Options
S0=100      #Stock price at time 0
K=80      #Strike (or exercise) price of the stock 
r=0.1      #risk free interest rate
sigma=0.25 #volatility (std deviation) of the underlying stock price
tau=1        #time (usually I put 1 year)
steps=12   #step size, i.e into how many steps you want your approximation discretized 

###OPTION 1 (starts)
#Package AmericanCallOpt

require(AmericanCallOpt)
call_price_am_bin<-am_call_bin(S0, K, r, sigma, tau, steps)
call_price_am_bin
###OPTION 1 (ends)

###OPTION 2 (starts)
#Matlab code translated into R - (better than using a package - see option 1)

AmericanCallDiv <- function(S0,K,r,tau,sigma,D=0,tauD=0,steps){
  
  f7 = 1
  dt = tau / steps
  v = exp(-r*dt )
  u = exp( sigma * sqrt ( dt ))
  d = 1 / u
  p = ( exp ( r * dt ) - d) / ( u - d )
  
  # adjust spot for dividend
  S0 = S0 - D * exp (- r * tauD )
  S = matrix(0, steps + 1 ,1)
  S[f7+0,] = S0 * d ^ steps
  
  for (j in 1:steps){
    S[f7 + j,] = S[f7 + j - 1,] * u / d   
  }
  
  # initialise option values at maturity ( period M )
  
  C = apply(cbind(S - K , 0),1,max)
  
  # step back through the tree
  
  for (i in (steps-1):0){
    for (j in 0:i){
      C[f7 +j] = v * ( p * C[f7 + j + 1] + (1 - p) * C[f7 + j] )
      S [f7 +j] = S[f7 + j] / d
      t = tau * i / steps 
      if (t > tauD){
        C[f7 + j] = max(C[f7 + j] , S[f7 + j,] - K )
      }
      else{
        C[f7 + j] = max(C[f7 + j] , S[f7 + j,] + D*exp(-r*(tauD-t)) - K) 
      }
    }
  }
  
  return(C0 = C[f7+0])
}
AmericanCallDiv(S0,K,r,tau,sigma,D=0,tauD=0,steps)

###OPTION 2 (ends)

#Note: Option 1 and 2 give exactly the same answer, which is great

###OPTION 3 (starts) 
#European call option just so you see that there exists a price difference
#between the two option styles

EU_blackscholes <-function(S0,K,r,sigma,tau,flag="c"){
  # s0 . . . i n i t i a l a s s e t p r i c e
  # K . . . s t r i k e P r ice
  # r . . . r i s k −f r e e i n t e r e s t r a t e
  # s igma . . . v o l a t i l i t y
  # T . . . t ime h o r i z o n
  # f l a g . . . c f o r c a l l o p t i o n s
  d1 <- (log(S0/K)+(r+sigma^2/2 )*tau)/(sigma*sqrt(tau))
  d2 <- d1-sigma*sqrt(tau)
  if (flag=="c"){
    S0*pnorm(d1)-K*exp(-r*tau)*pnorm(d2)
  } else{
    K*exp(-r*tau)*pnorm(-d2)-S0*pnorm(-d1)
  }
}
EU_blackscholes(S0,K,r,sigma,tau)

###OPTION 3 (ends)

#Now that seems quite deterministic, let's try some other model:

###OPTION 4 (starts)
#A.2. R Codes for the Implementation of the LSM Algorithm

n=1000 #I suppose the number of simulations?
d=1000 #used in th loop (number of iter?)

simul_AMERICAN_LSM <- function(n,d,S0,K,sigma,r,tau){
  S0<-S0/K
  dt<-tau/d
  z<-rnorm( n )
  s.t<-S0*exp((r - 1/2*sigma ^ 2)*tau+sigma*z*(tau^ 0.5))
  s.t[(n + 1):(2*n)]<-S0*exp((r - 1/2*sigma ^ 2 )*tau-sigma*z*(tau^ 0.5))
  CC<-pmax(1-s.t,0)

  payoffeu<-exp(-r*T)*(CC[1:n]+CC[(n+1):(2*n)])/2*K
  euprice<-mean(payoffeu)
  
  for (k in (d-1):1){
    z<-rnorm(n)
    mean<-(log(S0)+k*log(s.t[1 : n] ) ) /( k+1)
    vol<-( k*dt/( k + 1 ) ) ^ 0.5*z
    s.t_1<-exp(mean+sigma*vol)
    mean<-(log(S0)+k*log(s.t[(n + 1): ( 2*n )]) ) /(k+1)
    s.t_1[(n + 1): (2*n)]<-exp(mean-sigma*vol)
    CE<-pmax(1-s.t_1,0)
    idx <- (1 : ( 2 *n ) ) [CE>0]
    discountedCC<- CC[ idx ] *exp(-r*dt)
    basis1<-exp(-s.t_1[idx] /2 )
    basis2<-basis1*(1-s.t_1[ idx ] )
    basis3<-basis1*(1-2*s.t_1[ idx ]+(s.t_1[ idx ] ^ 2 ) /2 )
    p<-glm(discountedCC~basis1+basis2+basis3)$coefficients
    estimatedCC<-p[1]+ p[2] * basis1+p[ 3 ] * basis2+p[ 4 ] * basis3
    EF<-rep( 0 , 2*n)
    EF[ idx ]<-(CE[ idx ]>estimatedCC)
    CC<-(EF==0)*CC*exp(-r*dt)+(EF==1)*CE
    s.t<-s.t_1
  }
  
  payoff<-exp(-r*dt)*(CC[ 1 : n]+CC[ ( n + 1 ): ( 2*n ) ] ) /2
  usprice<-mean(payoff*K)
  error<-1.96*sd(payoff*K)/sqrt(n)
  earlyex<-usprice - euprice
  data.frame(usprice, error, euprice) #,premium? (add it if you can find the column)
}
simul_AMERICAN_LSM(n,d,S0,K,sigma,r,tau)  

#definitely a big price difference here!

###OPTION 3 (ends)


###OPTION 4 (starts)
#American Put Options

AmericanPutDiv <- function(S0,K,r,tau,sigma,D=0,tauD=0,steps){
  
  f7 = 1
  dt = tau / steps
  v = exp(-r*dt )
  u = exp( sigma * sqrt ( dt ))
  d = 1 / u
  p = ( exp ( r * dt ) - d) / ( u - d )
  
  # adjust spot for dividend
  S0 = S0 - D * exp (- r * tauD )
  S = matrix(0, steps + 1 ,1)
  S[f7+0,] = S0 * d ^ steps
  
  for (j in 1:steps){
    S[f7 + j,] = S[f7 + j - 1,] * u / d   
  }
  
  # initialise option values at maturity ( period M )
  
  C = apply(cbind(K - S , 0),1,max)
  
  # step back through the tree
  
  for (i in (steps-1):0){
    for (j in 0:i){
      C[f7 +j] = v * ( p * C[f7 + j + 1] + (1 - p) * C[f7 + j] )
      S [f7 +j] = S[f7 + j] / d
      t = tau * i / steps 
      if (t > tauD){
        C[f7 + j] = max(C[f7 + j] , S[f7 + j,] - K )
      }
      else{
        C[f7 + j] = max(C[f7 + j] , S[f7 + j,] + D*exp(-r*(tauD-t)) - K) 
      }
    }
  }
  
  return(C0 = C[f7+0])
}
AmericanPutDiv(S0,K,r,tau,sigma,D=0,tauD=0,steps)

###Option 4 (ends)

###Option 5 (starts)

install.packages("fOptions")
library(fOptions)

###Option 5 (ends)

###PLOTS AND CO (starts)
#A.3. R Codes for the Exact Early Exercise Boundary of an American Put
#at time td−1 by Newton-Raphson Method
newton_raphson<-function (S0,K,r,sigma,tau,steps){
  s<-S0
  dt <- tau/steps
  repeat{
    d1=(log ( s/K)+( r+sigma ^2/2 )*dt )/( sigma*sqrt(dt)) 
    d2=d1-sigma*sqrt(dt)
    f<-(K*exp(-r*dt)*pnorm(-d2)-s*pnorm(-d1)-(K-s))
    fderivative<-((s*sigma*sqrt(dt)) ^ (-1) )*(-K*exp(-r*dt)*dnorm(d2)+s*dnorm(d1))+pnorm(d1)
    new.s<-s-f/fderivative
    conv<-abs(new.s-s)
    if((conv/abs(s))<1e-8) break
    s<-new.s
  }
  s
}
newton_raphson(S0,K,r,sigma,tau,steps)


###Yerik's code

###
sim.Price_Bin <- function(N=10, u=3/2, d=2/3, x0 =1){
  
  x <- matrix(0, nrow = N+1, ncol = N+1)
  
  for(i in 1: (N+1)){
    
    for(j in 1:i){
      
      x[i, j]=x0*u^(j-1)*d^(i-j)
    }
    
  }
  return(x)
}
####

####
sim.Prob_Bin <- function(N=10, p = 0.4121){
  
  pi <- matrix(0, nrow = N+1, ncol = N+1)
  
  for(i in 1: (N+1)){
    
    for(j in 1:i){
      
      pi[i,j] <- dbinom(x = j-1, size = i-1, prob = p )
    }
    
  }
  
  return(pi)
}
####


g <- function(x, strike = 1){
  
  return(pmax(0, strike - x))
  
}


P <- function(state, p = 0.4121){
  return(sum(state*c(1-p, p))) 
}


N <- 3
x_strike <- 1
alpha <- 0.99
p = 0.4121
x0=1

x <- sim.Price_Bin(N=N, x0 = x0)
pi <- sim.Prob_Bin(N=N, p = p)



J <- matrix(0, nrow = N+1, ncol = N+1)


J[N+1,] <- g(x=x[N+1,])

for(i in N:1){
  
  for(j in 1:i){
    J[i,j] <- pmax(g(x=x[i,j]), alpha * P(J[i+1,c(j,j+1)], p))
    
  }
  
}


##Q value


Q <- matrix(0, nrow = N+1, ncol = N+1)

Q[N+1,] <- g(x[N+1,])

for(j in 1:N){
  
  Q[N, j] <- alpha * P(Q[N+1,c(j,j+1)], p)
  
}

for(i in (N-1):1){
  for(j in 1:i){
    Q[i,j] <- alpha * P(pmax(g(x[i+1,c(j, j+1)]), Q[i+1,c(j,j+1)]), p)
    
  }
  
}

J_Q <- pmax(g(x[1,1]), Q[1,1])


##Regression with fictitious pi

dim <- 3

phi1 <- function(x){
  
  return(rep(1, length(x)))
}

phi2 <- function(x){
  
  return(x)
}

phi3 <- function(x){
  
  return(x^2)
}


phi <- function(x, dim = 3, base = list(phi1, phi2, phi3)){
  
  base_Mat <- matrix(0, nrow = length(x), ncol = dim)
  for( j in 1:dim){
    base_Mat[,j] <- base[[j]](x)
    
  }
  return(base_Mat)
}

r <- matrix(0, ncol = N, nrow = dim )


S <- seq(from = 0.1, to = 2, by = 0.1)

u = 3/2
d = 2/3

Q1 <- matrix(0, nrow = N, ncol = length(S))

for(j in 1:length(S)){
  
  Q1[N, j] <- alpha * P(g(c(S[j]*d, S[j]*u)),p)
  
}

r[, N] <- solve((t(phi(S)) %*% phi(S)), t(phi(S)) %*% Q1[N,])


for (i in (N-1):1){
  
  for(j in 1:length(S)){
    
    
    Q1[i, j] <- alpha*P(pmax(g(c(S[j]*d, S[j]*u)), 
                             c(phi(S[j]*d) %*% r[,i+1], phi(S[j]*u) %*% r[,i+1])))
    
  }
  
  r[, i] <- solve((t(phi(S)) %*% phi(S)), t(phi(S)) %*% Q1[i,])
  
}

Q0_1 <- phi(S[10]) %*% r[,1]
##regression real





