rm(list = ls())

# maximum x for truncation
max_x <- 60

# ZTGID pmf
pmf_ztgid <- function(x, theta){
  
  const <- (theta^3/(theta^3 + 2)) /
    (1 - (theta^3/(theta^3 + 2)) *
       (theta/(theta + 1) + 2/(theta + 1)^3))
  
  k <- 0:x
  
  inner_sum <- sum((-1)^k * choose(x, k) *
                     (theta/(theta + k + 1) +
                        2/(theta + k + 1)^3))
  
  return(const * inner_sum)
}

# MGF function
mgf_ztgid <- function(t, theta){
  
  x_vals <- 1:max_x
  pmf_vals <- sapply(x_vals, pmf_ztgid, theta = theta)
  
  sum(exp(t * x_vals) * pmf_vals)
}

# t values
t_vals <- seq(-0.5, 0.5, length = 200)

# theta values
theta_vals <- c(1,1.5,2,2.5)

# compute mgf values
mgf_matrix <- sapply(theta_vals, function(th){
  sapply(t_vals, mgf_ztgid, theta = th)
})

# plot
plot(t_vals, mgf_matrix[,1], type="l", lwd=2,
     xlab="t", ylab="MGF",
     main="MGF of ZTGID for different theta")

lines(t_vals, mgf_matrix[,2], col=2, lwd=2)
lines(t_vals, mgf_matrix[,3], col=3, lwd=2)
lines(t_vals, mgf_matrix[,4], col=4, lwd=2)

legend("topleft",
       legend=c("theta = 1","theta = 1.5","theta = 2","theta = 2.5"),
       col=1:4, lwd=2)
