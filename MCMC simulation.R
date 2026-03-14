#n=100
rm(list = ls())

theta_values <- c(1,1.5,2,2.5,3,3.5,4)

n <- 100
B <- 1000
max_x <- 50


# PMF of ZTGID
pmf_ztgid <- function(x, theta){
  
  const <- theta^3/(theta^3 + 2)
  
  s <- 0
  for(k in 0:x){
    s <- s + (-1)^k * choose(x,k) *
      ( theta/(theta+k+1) + 2/(theta+k+1)^3 )
  }
  
  num <- const * s
  den <- 1 - const * ( theta/(theta+1) + 2/(theta+1)^3 )
  
  return(num/den)
}


# Score function (MLE equation)
score_fun <- function(theta,x){
  
  #n <- length(x)
  total <- 0
  
  for(i in 1:n){
    
    xi <- x[i]
    
    part1 <- 0
    part2 <- 0
    
    for(k in 0:xi){
      
      part1 <- part1 + (-1)^k*choose(xi,k)*
        ( theta/(theta+k+1) + 2/(theta+k+1)^3 )
      
      part2 <- part2 + (-1)^k*choose(xi,k)*
        ( (k+1)/(theta+k+1)^2 - 6/(theta+k+1)^4 )
      
    }
    
    total <- total + (6/(theta*(theta^3+2)))*part1 + part2
  }
  
  return(total)
}


# Storage for results
results <- data.frame(True_theta = theta_values,
                      Mean_theta = NA,
                      SE_theta = NA)


# Show multiple histograms
par(mfrow=c(3,3))


for(t in 1:length(theta_values)){
  
  theta_true <- theta_values[t]
  
  theta_est <- numeric(B)
  
  
  for(b in 1:B){
    
    # Possible x values
    x_vals <- 1:max_x
    
    # Probabilities
    p_vals <- sapply(x_vals, pmf_ztgid, theta=theta_true)
    
    
    
    # Generate MCMC 
    x <- numeric(n)
    x[1] <- sample(x_vals,1)
    
    for(i in 2:n){
      
      current <- x[i-1]
      proposal <- sample(x_vals,1)
      
      alpha <- min(1, p_vals[proposal]/p_vals[current])
      
      if(runif(1) < alpha){
        x[i] <- proposal
      }else{
        x[i] <- current
      }
      
    }
    
    
    # Solve MLE using bisection (uniroot)
    root <- tryCatch(
      uniroot(function(th) score_fun(th,x), c(0.5,10))$root,
      error=function(e) NA
    )
    
    theta_est[b] <- root
    
  }
  
  
  # Mean and Standard Error
  results$Mean_theta[t] <- mean(theta_est, na.rm=TRUE)
  results$SE_theta[t] <- sd(theta_est, na.rm=TRUE)
  
  
  # Histogram
  hist(theta_est,
       main=paste("Histogram of Theta Estimates (True θ =",theta_true,")"),
       xlab="Theta estimates",
       col="skyblue",
       border="black")
  
}


# Print results
results


#n=50
rm(list = ls())

theta_values <- c(1,1.5,2,2.5,3,3.5,4)

n <- 50
B <- 1000
max_x <- 50


# PMF of ZTGID
pmf_ztgid <- function(x, theta){
  
  const <- theta^3/(theta^3 + 2)
  
  s <- 0
  for(k in 0:x){
    s <- s + (-1)^k * choose(x,k) *
      ( theta/(theta+k+1) + 2/(theta+k+1)^3 )
  }
  
  num <- const * s
  den <- 1 - const * ( theta/(theta+1) + 2/(theta+1)^3 )
  
  return(num/den)
}


# Score function (MLE equation)
score_fun <- function(theta,x){
  
  #n <- length(x)
  total <- 0
  
  for(i in 1:n){
    
    xi <- x[i]
    
    part1 <- 0
    part2 <- 0
    
    for(k in 0:xi){
      
      part1 <- part1 + (-1)^k*choose(xi,k)*
        ( theta/(theta+k+1) + 2/(theta+k+1)^3 )
      
      part2 <- part2 + (-1)^k*choose(xi,k)*
        ( (k+1)/(theta+k+1)^2 - 6/(theta+k+1)^4 )
      
    }
    
    total <- total + (6/(theta*(theta^3+2)))*part1 + part2
  }
  
  return(total)
}


# Storage for results
results <- data.frame(True_theta = theta_values,
                      Mean_theta = NA,
                      SE_theta = NA)


# Show multiple histograms
par(mfrow=c(3,3))


for(t in 1:length(theta_values)){
  
  theta_true <- theta_values[t]
  
  theta_est <- numeric(B)
  
  
  for(b in 1:B){
    
    # Possible x values
    x_vals <- 1:max_x
    
    # Probabilities
    p_vals <- sapply(x_vals, pmf_ztgid, theta=theta_true)
    
    
    
    # Generate MCMC sample
    x <- numeric(n)
    x[1] <- sample(x_vals,1)
    
    for(i in 2:n){
      
      current <- x[i-1]
      proposal <- sample(x_vals,1)
      
      alpha <- min(1, p_vals[proposal]/p_vals[current])
      
      if(runif(1) < alpha){
        x[i] <- proposal
      }else{
        x[i] <- current
      }
      
    }
    
    
    # Solve MLE using bisection (uniroot)
    root <- tryCatch(
      uniroot(function(th) score_fun(th,x), c(0.5,10))$root,
      error=function(e) NA
    )
    
    theta_est[b] <- root
    
  }
  
  
  # Mean and Standard Error
  results$Mean_theta[t] <- mean(theta_est, na.rm=TRUE)
  results$SE_theta[t] <- sd(theta_est, na.rm=TRUE)
  
  
  # Histogram
  hist(theta_est,
       main=paste("Histogram of Theta Estimates (True θ =",theta_true,")"),
       xlab="Theta estimates",
       col="skyblue",
       border="black")
  
}


# Print results
results



#n=200
rm(list = ls())

theta_values <- c(1,1.5,2,2.5,3,3.5,4)

n <- 200
B <- 1000
max_x <- 50


# PMF of ZTGID
pmf_ztgid <- function(x, theta){
  
  const <- theta^3/(theta^3 + 2)
  
  s <- 0
  for(k in 0:x){
    s <- s + (-1)^k * choose(x,k) *
      ( theta/(theta+k+1) + 2/(theta+k+1)^3 )
  }
  
  num <- const * s
  den <- 1 - const * ( theta/(theta+1) + 2/(theta+1)^3 )
  
  return(num/den)
}


# Score function (MLE equation)
score_fun <- function(theta,x){
  
  #n <- length(x)
  total <- 0
  
  for(i in 1:n){
    
    xi <- x[i]
    
    part1 <- 0
    part2 <- 0
    
    for(k in 0:xi){
      
      part1 <- part1 + (-1)^k*choose(xi,k)*
        ( theta/(theta+k+1) + 2/(theta+k+1)^3 )
      
      part2 <- part2 + (-1)^k*choose(xi,k)*
        ( (k+1)/(theta+k+1)^2 - 6/(theta+k+1)^4 )
      
    }
    
    total <- total + (6/(theta*(theta^3+2)))*part1 + part2
  }
  
  return(total)
}


# Storage for results
results <- data.frame(True_theta = theta_values,
                      Mean_theta = NA,
                      SE_theta = NA)


# Show multiple histograms
par(mfrow=c(3,3))


for(t in 1:length(theta_values)){
  
  theta_true <- theta_values[t]
  
  theta_est <- numeric(B)
  
  
  for(b in 1:B){
    
    # Possible x values
    x_vals <- 1:max_x
    
    # Probabilities
    p_vals <- sapply(x_vals, pmf_ztgid, theta=theta_true)
    
    
    
    # Generate MCMC sample
    x <- numeric(n)
    x[1] <- sample(x_vals,1)
    
    for(i in 2:n){
      
      current <- x[i-1]
      proposal <- sample(x_vals,1)
      
      alpha <- min(1, p_vals[proposal]/p_vals[current])
      
      if(runif(1) < alpha){
        x[i] <- proposal
      }else{
        x[i] <- current
      }
      
    }
    
    
    # Solve MLE using bisection (uniroot)
    root <- tryCatch(
      uniroot(function(th) score_fun(th,x), c(0.5,10))$root,
      error=function(e) NA
    )
    
    theta_est[b] <- root
    
  }
  
  
  # Mean and Standard Error
  results$Mean_theta[t] <- mean(theta_est, na.rm=TRUE)
  results$SE_theta[t] <- sd(theta_est, na.rm=TRUE)
  
  
  # Histogram
  hist(theta_est,
       main=paste("Histogram of Theta Estimates (True θ =",theta_true,")"),
       xlab="Theta estimates",
       col="skyblue",
       border="black")
  
}


# Print results
results


#n=50
#MOME 

rm(list = ls())

theta_values <- c(1,1.5,2,2.5,3,3.5,4)

n <- 50
B <- 1000
max_x <- 50


# PMF of ZTGID
pmf_ztgid <- function(x, theta){
  
  const <- theta^3/(theta^3 + 2)
  
  s <- 0
  for(k in 0:x){
    s <- s + (-1)^k * choose(x,k) *
      ( theta/(theta+k+1) + 2/(theta+k+1)^3 )
  }
  
  num <- const * s
  den <- 1 - const * ( theta/(theta+1) + 2/(theta+1)^3 )
  
  return(num/den)
}


# MOME equation
mome_fun <- function(theta, xbar){
  
  (xbar-1)*theta^6 -
    (3*xbar-2)*theta^5 +
    (3*xbar-1)*theta^4 -
    (xbar+1)*theta^3 -
    6*xbar*theta^2 -
    2*xbar
  
}


# Result storage
results <- data.frame(True_theta = theta_values,
                      Mean_theta = NA,
                      SE_theta = NA)


# Multiple histograms
par(mfrow=c(3,3))


for(t in 1:length(theta_values)){
  
  theta_true <- theta_values[t]
  
  theta_est <- numeric(B)
  
  
  for(b in 1:B){
    
    # Possible values
    x_vals <- 1:max_x
    
    # Probabilities
    p_vals <- sapply(x_vals, pmf_ztgid, theta=theta_true)
    
    
    
    # MCMC sample generation
    x <- numeric(n)
    x[1] <- sample(x_vals,1)
    
    for(i in 2:n){
      
      current <- x[i-1]
      proposal <- sample(x_vals,1)
      
      alpha <- min(1, p_vals[proposal]/p_vals[current])
      
      if(runif(1) < alpha){
        x[i] <- proposal
      }else{
        x[i] <- current
      }
      
    }
    
    
    # Sample mean
    xbar <- mean(x)
    
    
    # Solve MOME equation
    root <- tryCatch(
      uniroot(function(th) mome_fun(th,xbar), c(0.1,10))$root,
      error=function(e) NA
    )
    
    
    theta_est[b] <- root
    
  }
  
  
  # Mean and SE
  results$Mean_theta[t] <- mean(theta_est, na.rm=TRUE)
  results$SE_theta[t] <- sd(theta_est, na.rm=TRUE)
  
  
  # Histogram
  hist(theta_est,
       main=paste("MOME Estimates (True θ =",theta_true,")"),
       xlab="Theta estimates",
       col="lightgreen",
       border="black")
  
}


# Final simulation table
results



#n=200
#MOME 

rm(list = ls())

theta_values <- c(1,1.5,2,2.5,3,3.5,4)

n <- 200
B <- 1000
max_x <- 50


# PMF of ZTGID
pmf_ztgid <- function(x, theta){
  
  const <- theta^3/(theta^3 + 2)
  
  s <- 0
  for(k in 0:x){
    s <- s + (-1)^k * choose(x,k) *
      ( theta/(theta+k+1) + 2/(theta+k+1)^3 )
  }
  
  num <- const * s
  den <- 1 - const * ( theta/(theta+1) + 2/(theta+1)^3 )
  
  return(num/den)
}


# MOME equation
mome_fun <- function(theta, xbar){
  
  (xbar-1)*theta^6 -
    (3*xbar-2)*theta^5 +
    (3*xbar-1)*theta^4 -
    (xbar+1)*theta^3 -
    6*xbar*theta^2 -
    2*xbar
  
}


# Result storage
results <- data.frame(True_theta = theta_values,
                      Mean_theta = NA,
                      SE_theta = NA)


# Multiple histograms
par(mfrow=c(3,3))


for(t in 1:length(theta_values)){
  
  theta_true <- theta_values[t]
  
  theta_est <- numeric(B)
  
  
  for(b in 1:B){
    
    # Possible values
    x_vals <- 1:max_x
    
    # Probabilities
    p_vals <- sapply(x_vals, pmf_ztgid, theta=theta_true)
    
    
    
    # MCMC sample generation
    x <- numeric(n)
    x[1] <- sample(x_vals,1)
    
    for(i in 2:n){
      
      current <- x[i-1]
      proposal <- sample(x_vals,1)
      
      alpha <- min(1, p_vals[proposal]/p_vals[current])
      
      if(runif(1) < alpha){
        x[i] <- proposal
      }else{
        x[i] <- current
      }
      
    }
    
    
    # Sample mean
    xbar <- mean(x)
    
    
    # Solve MOME equation
    root <- tryCatch(
      uniroot(function(th) mome_fun(th,xbar), c(0.1,10))$root,
      error=function(e) NA
    )
    
    
    theta_est[b] <- root
    
  }
  
  
  # Mean and SE
  results$Mean_theta[t] <- mean(theta_est, na.rm=TRUE)
  results$SE_theta[t] <- sd(theta_est, na.rm=TRUE)
  
  
  # Histogram
  hist(theta_est,
       main=paste("MOME Estimates (True θ =",theta_true,")"),
       xlab="Theta estimates",
       col="lightgreen",
       border="black")
  
}


# Final simulation table
results


