# PMF of parent GID

pmf_gid <- function(x, theta){
  
  const <- theta^3/(theta^3 + 2)
  
  s <- 0
  
  for(k in 0:x){
    s <- s + ((-1)^k) * choose(x,k) *
      ( theta/(theta+k+1) + 2/(theta+k+1)^3 )
  }
  
  return(const * s)
}

#Zero-truncated PMF

pmf_ztgid <- function(x, theta){
  
  p0 <- pmf_gid(0, theta)
  
  return( pmf_gid(x, theta) / (1 - p0) )
}

#CDF

cdf_ztgid <- function(x, theta){
  
  sum(sapply(1:x, function(t) pmf_ztgid(t, theta)))
}

# Generate values

x <- 1:20   # for adjust only

theta_vals <- c(1,1.5,2,2.5)

cdf_matrix <- sapply(theta_vals, function(th){
  sapply(x, function(xx) cdf_ztgid(xx, th))
})

# Plot

plot(x, cdf_matrix[,1],
     type="o",
     lwd=2,
     ylim=c(0,1),
     xlab="x",
     ylab = expression(paste("F(x,", theta, ")")),
     main="CDF of ZTGID")

lines(x, cdf_matrix[,2], type="o", lwd=2, col=2)
lines(x, cdf_matrix[,3], type="o", lwd=2, col=3)
lines(x, cdf_matrix[,4], type="o", lwd=2, col=4)

legend("bottomright",
       legend=c(expression(theta==1),
                expression(theta==1.5),
                expression(theta==2),
                expression(theta==2.5)),
       col=1:4,
       lwd=2)


