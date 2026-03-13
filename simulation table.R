rm(list=ls())

set.seed(123)

B <- 200
theta_vals <- c(0.1,0.4,1,1.5,3,4)
n_vals <- c(100,50,30)

pmf_ztgid <- function(x,theta){
  
  const <- theta^3/(theta^3+2)
  
  num <- sapply(x,function(xx){
    
    sum((-1)^(0:xx)*choose(xx,0:xx)*
          ( theta/(theta+(0:xx)+1) + 2/(theta+(0:xx)+1)^3 ))
  })
  
  p0 <- const*(theta/(theta+1)+2/(theta+1)^3)
  
  p <- const*num/(1-p0)
  
  p[p<0] <- 0
  p[!is.finite(p)] <- 0
  
  p/sum(p)
}

simulate_sample <- function(theta,n){
  
  x_vals <- 1:50
  probs <- pmf_ztgid(x_vals,theta)
  
  sample(x_vals,n,replace=TRUE,prob=probs)
}

MOME_est <- function(x){
  
  xbar <- mean(x)
  
  g <- function(theta){
    
    (xbar-1)*theta^6 -
      (3*xbar-2)*theta^5 +
      (3*xbar-1)*theta^4 -
      (xbar+1)*theta^3 -
      6*xbar*theta^2 -
      2*xbar
  }
  
  g_prime <- function(theta){
    
    6*(xbar-1)*theta^5 -
      5*(3*xbar-2)*theta^4 +
      4*(3*xbar-1)*theta^3 -
      3*(xbar+1)*theta^2 -
      12*xbar*theta
  }
  
  theta <- 2
  
  for(i in 1:100){
    
    theta_new <- theta-g(theta)/g_prime(theta)
    
    if(abs(theta_new-theta)<1e-8) break
    
    theta <- theta_new
  }
  
  theta
}

MLE_est <- function(x){
  
  n <- length(x)
  
  loglik <- function(theta){
    
    if(theta<=1) return(-1e10)
    
    const <- theta^3/(theta^3+2)
    
    p0 <- const*(theta/(theta+1)+2/(theta+1)^3)
    
    part1 <- n*log(const)
    part3 <- -n*log(1-p0)
    
    part2 <- 0
    
    for(i in 1:n){
      
      xi <- x[i]
      
      inner <- sum((-1)^(0:xi)*choose(xi,0:xi)*
                     ( theta/(theta+(0:xi)+1)+2/(theta+(0:xi)+1)^3 ))
      
      if(inner<=0) return(-1e10)
      
      part2 <- part2+log(inner)
    }
    
    part1+part2+part3
  }
  
  optimize(function(t)-loglik(t),c(1.01,20))$minimum
}

results <- data.frame()

for(theta in theta_vals){
  
  row <- c(theta)
  
  for(n in n_vals){
    
    mome_vec <- numeric(B)
    mle_vec <- numeric(B)
    
    for(b in 1:B){
      
      x <- simulate_sample(theta,n)
      
      mome_vec[b] <- tryCatch(MOME_est(x),error=function(e) NA)
      
      mle_vec[b] <- tryCatch(MLE_est(x),error=function(e) NA)
    }
    
    mome_mean <- mean(mome_vec,na.rm=TRUE)
    mome_se <- sd(mome_vec,na.rm=TRUE)
    
    mle_mean <- mean(mle_vec,na.rm=TRUE)
    mle_se <- sd(mle_vec,na.rm=TRUE)
    
    row <- c(row,mome_mean,mome_se,mle_mean,mle_se)
  }
  
  results <- rbind(results,row)
  
  cat("Completed theta =",theta,"\n")
}

colnames(results) <- c(
  "Theta",
  "MOME100","SE100","MLE100","SE100_MLE",
  "MOME50","SE50","MLE50","SE50_MLE",
  "MOME30","SE30","MLE30","SE30_MLE"
)

print(results)

