
theta_true <- 3
n <- 100

pmf_ztgid <- function(x, theta){
  
  const <- theta^3/(theta^3 + 2)
  
  numerator <- 0
  for(k in 0:x){
    numerator <- numerator +
      (-1)^k * choose(x,k) *
      ( theta/(theta+k+1) +
          2/(theta+k+1)^3 )
  }
  
  p0_0 <- const * ( theta/(theta+1) +
                      2/(theta+1)^3 )
  
  val <- const * numerator / (1 - p0_0)
  return(val)
}

x_vals <- 1:50
probs <- sapply(x_vals, pmf_ztgid, theta = theta_true)

probs[probs < 0] <- 0
probs[!is.finite(probs)] <- 0
probs <- probs / sum(probs)

set.seed(123)
x_sample <- sample(x_vals,
                   size = n,
                   replace = TRUE,
                   prob = probs)

xbar <- mean(x_sample)

g <- function(theta){
  (xbar - 1)*theta^6 -
    (3*xbar - 2)*theta^5 +
    (3*xbar - 1)*theta^4 -
    (xbar + 1)*theta^3 -
    6*xbar*theta^2 -
    2*xbar
}

g_prime <- function(theta){
  6*(xbar - 1)*theta^5 -
    5*(3*xbar - 2)*theta^4 +
    4*(3*xbar - 1)*theta^3 -
    3*(xbar + 1)*theta^2 -
    12*xbar*theta
}

newton_mome <- function(theta_init,
                        tol = 1e-10,
                        max_iter = 100){
  
  theta <- theta_init
  
  for(i in 1:max_iter){
    
    deriv <- g_prime(theta)
    if(abs(deriv) < 1e-8)
      stop("Derivative too small in MOME NR")
    
    theta_new <- theta - g(theta)/deriv
    
    if(abs(theta_new - theta) < tol)
      return(theta_new)
    
    theta <- theta_new
  }
  
  stop("MOME NR did not converge")
}

theta_MOME <- newton_mome(2)


loglik <- function(theta){
  
  if(!is.finite(theta) || theta <= 1)
    return(-1e12)
  
  const <- theta^3/(theta^3 + 2)
  part1 <- n * log(const)
  
  p0_0 <- const * ( theta/(theta+1) +
                      2/(theta+1)^3 )
  
  denom <- 1 - p0_0
  if(denom <= 0) return(-1e12)
  
  part3 <- -n * log(denom)
  
  part2 <- 0
  
  for(i in 1:n){
    
    xi <- x_sample[i]
    inner <- 0
    
    for(k in 0:xi){
      inner <- inner +
        (-1)^k * choose(xi,k) *
        ( theta/(theta+k+1) +
            2/(theta+k+1)^3 )
    }
    
    if(inner <= 0) return(-1e12)
    
    part2 <- part2 + log(inner)
  }
  
  return(part1 + part2 + part3)
}


score <- function(theta, h = 1e-4){
  if(theta <= 1) return(0)
  (loglik(theta + h) - loglik(theta - h)) / (2*h)
}



fixed_point_mle <- function(theta_init,
                            alpha = 1e-5,
                            tol = 1e-6,
                            max_iter = 500,
                            upper_bound = 20){
  
  theta <- theta_init
  
  for(i in 1:max_iter){
    
    s <- score(theta)
    theta_new <- theta + alpha * s
    
    if(theta_new <= 1) theta_new <- 1.05
    if(theta_new > upper_bound) theta_new <- upper_bound
    
    if(loglik(theta_new) < loglik(theta))
      theta_new <- theta
    
    if(abs(theta_new - theta) < tol)
      return(theta_new)
    
    theta <- theta_new
  }
  
  return(theta)
}

theta_MLE_FP <- fixed_point_mle(2)

cat("Sample size (n):", n, "\n")
cat("True theta      :", theta_true, "\n")
cat("MOME (NR)       :", theta_MOME, "\n")
cat("MLE (Fixed Point):", theta_MLE_FP, "\n")


theta_true <- 3


theta_MOME <- 2.815894
theta_MLE  <- 1.937014


rel_MOME <- abs(theta_MOME - theta_true) / theta_true * 100
rel_MLE  <- abs(theta_MLE  - theta_true) / theta_true * 100

cat("Relative difference (MOME):", round(rel_MOME, 4), "%\n")
cat("Relative difference (MLE) :", round(rel_MLE, 4), "%\n")
