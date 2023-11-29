#### e-BH #


K = 40 #### total number of experiments  
K0 = 30 ##### number of nulls 
mu = 2 #### signal strength
delta = 4  #### martingale construction  
alpha = c(0.1,0.05,0.02) #### e-BH FDR level 
penalty = sum(1/(1:K)) ##### penalty paid by BY
alphaBY=alpha/penalty #### adjusted threshold for BY
rho = 0 #### correlation, non-negative
sqr = sqrt(rho) 
setting = 1 ##### 1 for non-negative correlation, 2 for -1/K correaltion, 3 for banded -0.5 correlation
gam <- 0.01


### Define the number of rejections of e-BH or BH using input e-values
ReBH <- function(TE,alpha){
  p_vals_adjusted <- p.adjust(1 / TE, method = "BH")
  rejected_indices <- which(p_vals_adjusted <= alpha)
  return(length(rejected_indices))
}


Truedis <- function(TE,alpha,B){
  TEnB <- rbind(TE,B)
  TEnBsort <-  TEnB[,order(TEnB[1,],decreasing=TRUE)]
  r <- ReBH(TE,alpha)
  if (r == 0){return(0)}
  else{return(sum(TEnBsort[2,1:r]))}
}



#### START EXPERIMENTS 




N=10000  # repeat times!

start_time <- Sys.time() 
set.seed(5) 

#GG <- matrix(rep(0, len=15),nrow=5) ### record final results
HH <- matrix(rep(0, len=99),nrow=11) ### record final results



colnames(HH) <- c(alpha[1],"FDP%", "PWR", alpha[2],"FDP%", "PWR", alpha[3], "FDP%", "PWR")
rownames(HH) <- c("BH", "e-BH PRDS", "BY", "e-BH AD", "e-BH", "e-BH Story", "e-BH AD Storey", "e-BH PRDS Storey", "e-BH TS", "e-BH AD TS", "e-BH PRDS TS") 

pi_0_vec <- c()
pi_0_vec2 <- c()

for (h in 1:N)
{
  
  B = c(rep(0,len=K0),rep(1,len=K-K0))
  B <- sample(B) 
  
  
  ##### Data generation
  
  
  X<- rep(0, len=K)
  if (setting == 1){
    Z <- rnorm(1)
    Y <- rnorm(K) 
    X <- sqrt(1-sqr^2)*Y + sqr*Z + B*mu
  }
  else{
    if (setting == 2){
      Y <- rnorm(K)
      X <- sqrt(K/(K-1))*(-Y + sum(Y)/K)
      X <- X+ B*mu
    }
    else{
      Y <- rnorm(K+1)
      for (l in 1:K){
        X[l] <-sqrt(1/2) *( Y[l]-Y[l+1]) 
        X <- X+B*mu
      }
    }
  }
  
  
  P <- 1-pnorm(X)
  E <- exp(delta*X-delta^2/2)
  PE <-1/P
  
  
  #### SOLVING FOR BOOSTING FACTORS.
  
  d=delta
  #boost1<-rep(0,len=3)
  #boost2<-rep(0,len=3)
  H <- matrix(rep(0,len=99), nrow=11) ### data recorder
  
  
  for (i in 1:3) {
    a=alpha[i]
    p<-rep(0,len=K)
    q<-rep(0,len=K) 
    boostAD <- function(b){
      p[1]=1-pnorm(log(K/1/a/b)/d + d/2)
      q[1]=K*p[1]
      for (k in 2:K){ 
        p[k]= pnorm(log(K/(k-1)/a/b)/d + d/2) -pnorm(log(K/k/a/b)/d + d/2)
        q[k]=p[k]*K/k
      }
      return(sum(q)-a) 
    }
    boost1<-uniroot(boostAD, c(1,1000))$root  ##### AD boost
    
    boostPR <- function(b){
      t=K/seq(1:K)
      z=max(t*pnorm(-d/2-log(t)/d + log(a  * b )/d))-a
      return(z)
    }
    boost2 <- uniroot(boostPR, c(1,1000))$root ##### PRDS boost
    
    #####  boost1=1.68 ##### AD boost for delta=3 and alpha=0.05, K=200
    #####  boost2=7.88 ##### PRDS boost for delta=3 and alpha=0.05, K=200
    
    E_sorted <- sort(E)
    k_hat <- 0
    for (k in (1:length(E))) {
      first_i <- E_sorted[1:k]
      if (mean(first_i) <= 1+gam) {
        k_hat <- k
      } else {
        break
      }
    }
    pi_0_hat <- (1 + k_hat) / length(E)
    
    beta <- alpha[i] / 10
    N_e_vals <- length(E)
    k_hat2 <- length(E[1/E > beta / N_e_vals])
    pi_0_hat2 <- (k_hat2) / length(E)
    
    pi_0_vec[h] <- pi_0_hat
    pi_0_vec2[h] <- pi_0_hat2
    
    Eb1=boost1*E
    Eb2=boost2*E 
    
    R1 <- ReBH(E,alpha[i])
    R2 <- ReBH(PE,alpha[i])
    R3 <- ReBH(PE,alphaBY[i])
    R4 <- ReBH(Eb1,alpha[i])
    R5 <- ReBH(Eb2,alpha[i])
    R6 <- ReBH(E,alpha[i] / pi_0_hat)
    R7 <- ReBH(Eb1,alpha[i] / pi_0_hat)
    R8 <- ReBH(Eb2,alpha[i] / pi_0_hat)
    R9 <- ReBH(E,(alpha[i] - beta) / pi_0_hat2)
    R10 <- ReBH(Eb1,(alpha[i] - beta) / pi_0_hat2)
    R11 <- ReBH(Eb2,(alpha[i] - beta) / pi_0_hat2)
    
    TD1 <- Truedis(E,alpha[i],B)
    TD2 <- Truedis(PE,alpha[i],B)
    TD3 <- Truedis(PE,alphaBY[i],B)
    TD4 <- Truedis(Eb1,alpha[i],B)
    TD5 <- Truedis(Eb2,alpha[i],B)
    TD6 <- Truedis(E,alpha[i] / pi_0_hat,B)
    TD7 <- Truedis(Eb1,alpha[i] / pi_0_hat,B)
    TD8 <- Truedis(Eb2,alpha[i] / pi_0_hat,B)
    TD9 <- Truedis(E,(alpha[i] - beta) / pi_0_hat2,B)
    TD10 <- Truedis(Eb1,(alpha[i] - beta) / pi_0_hat2,B)
    TD11 <- Truedis(Eb2,(alpha[i] - beta) / pi_0_hat2,B)
    
    H[1,3*i-2]<-R2
    H[2,3*i-2]<-R5
    H[3,3*i-2]<-R3
    H[4,3*i-2]<-R4
    H[5,3*i-2]<-R1
    H[6,3*i-2]<-R6
    H[7,3*i-2]<-R7
    H[8,3*i-2]<-R8
    H[9,3*i-2]<-R9
    H[10,3*i-2]<-R10
    H[11,3*i-2]<-R11
    
    H[1,3*i-1]<-((R2-TD2)/max(R2,1))*100
    H[2,3*i-1]<-((R5-TD5)/max(R5,1))*100
    H[3,3*i-1]<-((R3-TD3)/max(R3,1))*100
    H[4,3*i-1]<-((R4-TD4)/max(R4,1))*100
    H[5,3*i-1]<-((R1-TD1)/max(R1,1))*100
    H[6,3*i-1]<-((R6-TD6)/max(R6,1))*100
    H[7,3*i-1]<-((R7-TD7)/max(R7,1))*100
    H[8,3*i-1]<-((R8-TD8)/max(R8,1))*100
    H[9,3*i-1]<-((R9-TD9)/max(R9,1))*100
    H[10,3*i-1]<-((R10-TD10)/max(R10,1))*100
    H[11,3*i-1]<-((R11-TD11)/max(R11,1))*100
    
    H[1,3*i]<-(TD2 / (max(K - K0,1)))*100
    H[2,3*i]<-(TD5 / (max(K - K0,1)))*100
    H[3,3*i]<-(TD3 / (max(K - K0,1)))*100
    H[4,3*i]<-(TD4 / (max(K - K0,1)))*100
    H[5,3*i]<-(TD1 / (max(K - K0,1)))*100
    H[6,3*i]<-(TD6 / (max(K - K0,1)))*100
    H[7,3*i]<-(TD7 / (max(K - K0,1)))*100
    H[8,3*i]<-(TD8 / (max(K - K0,1)))*100
    H[9,3*i]<-(TD9 / (max(K - K0,1)))*100
    H[10,3*i]<-(TD10 / (max(K - K0,1)))*100
    H[11,3*i]<-(TD11 / (max(K - K0,1)))*100
  }
  HH = HH+H
} 
HH<-t(round(t(HH/N),digit=c(1,2,1,2,1,2))) ## averaging 

#G=t(matrix(c( 
#  R2, TD2,  (1-TD2/max(R2,1))*100,
#  R5, TD5,   (1-TD5/max(R5,1))*100,
#  R3, TD3,  (1-TD3/max(R3,1))*100,
#  R4, TD4,   (1-TD4/max(R4,1))*100,
#  R1, TD1,   (1-TD1/max(R1,1))*100
#), nrow=3))
#colnames(G) <- c("#dis", "true","FDP%" )

#rownames(G) <- c("BH", "e-BH PRDS", "BY", "e-BH AD", "e-BH") 

#GG <- GG+G ## summing results  
#print(h)
#
#}
#GG<-t(round(t(GG/N),digit=c(1,2,2))) ## averaging 

Para<-t(matrix(t(c("K","K0", "rho", "N","mu", "setting", K, K0,rho,N,mu,setting)),nrow=6))

end_time <- Sys.time()
end_time - start_time  

HH
Para
