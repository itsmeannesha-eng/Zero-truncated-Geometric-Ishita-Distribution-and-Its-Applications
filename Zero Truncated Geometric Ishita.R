
rm(list = ls())
graphics.off()

par(mar=c(4,4,2,1))

p_geom_ishita <- function(x, theta) {
  
  out <- numeric(length(x))
  
  for (i in seq_along(x)) {
    xx <- x[i]
    s <- 0
    
    for (k in 0:(xx)) {
      s <- s + (-1)^k * choose(xx, k) *
        (theta / (theta + k + 1) + 2 / (theta + k + 1)^3)
    }
    
    out[i] <- (theta^3 / (theta^3 + 2)) * s
  }
  
  out
}


p_zt_geom_ishita <- function(x, theta) {
  
  p0 <- (theta^3 / (theta^3 + 2)) *
    (theta / (theta + 1) + 2 / (theta + 1)^3)
  
  p_geom_ishita(x, theta) / (1 - p0)
}

x <- 1:15
theta_vals <- c(1, 1.5, 2, 2.5)

par(mfrow = c(2,2))

for(theta in theta_vals){
  
  probs <- p_zt_geom_ishita(x, theta)
  
  plot(x, probs,
       type="h",
       lwd=3,
       xlab="x",
       ylab="P(x)",
       main = bquote(theta == .(theta)),
       ylim=c(0,max(probs)))
  
 points(x, probs, pch=19)
}

par(mfrow=c(1,1))


theta <- seq(4.2,8,by=0.05)

par(mfrow=c(2,2), mar=c(4,4,2,1))



mean_term <- (theta^3/(theta^3+2))*
  (theta/(theta-1) + 2/(theta-1)^3)

mu2 <- (2*theta^3/(theta^3+2))*
  (2/(theta-2)^3 + theta/(theta-2)) -
  (theta^6/(theta^3+2)^2)*
  (theta/(theta-1) + 2/(theta-1)^3)^2

variance <- mu2



CV <- sqrt(variance)/mean_term

plot(theta, CV,
     type="l",
     lwd=3,
     col="blue",
     xlab=expression(theta),
     ylab="CV",
     main="Coefficient of Variation")
grid()



ID <- variance/mean_term

plot(theta, ID,
     type="l",
     lwd=3,
     col="red",
     xlab=expression(theta),
     ylab="Index",
     main="Index of Dispersion")
grid()



mu3 <- (theta^3/(theta^3+2))*
  (
    6*(theta/(theta-3) + 2/(theta-3)^3)
    - 6*(theta/(theta-2) + 2/(theta-2)^3)
    + (theta/(theta-1) + 2/(theta-1)^3)
  )
-
  (4*theta^3/(theta^3+2))*
  (2/(theta-2)^3 + theta/(theta-2))
-
  (2*theta^6/(theta^3+2)^2)*
  (theta/(theta-1) + 2/(theta-1)^3)^3


skewness <- mu3/(mu2^(3/2))

plot(theta, skewness,
     type="l",
     lwd=3,
     col="darkgreen",
     xlab=expression(theta),
     ylab="Skewness",
     main="Skewness")
grid()



mu4 <- (theta^3/(theta^3+2))*
  (
    24*(theta/(theta-4) + 2/(theta-4)^3)
    - 36*(theta/(theta-3) + 2/(theta-3)^3)
    + 14*(theta/(theta-2) + 2/(theta-2)^3)
    - (theta/(theta-1) + 2/(theta-1)^3)
  )

- 4*(theta^3/(theta^3+2))*
  (
    6*(theta/(theta-3) + 2/(theta-3)^3)
    - 6*(theta/(theta-2) + 2/(theta-2)^3)
    + (theta/(theta-1) + 2/(theta-1)^3)
  )*mean_term

+ 6*((2*theta^3/(theta^3+2))*
       (2/(theta-2)^3 + theta/(theta-2)))^2

- 3*(theta^12/(theta^3+2)^4)*
  (theta/(theta-1) + 2/(theta-1)^3)^4


kurtosis <- mu4/(mu2^2)

plot(theta, kurtosis,
     type = "l",
     lwd = 3,
     col = "purple",
     xlab = expression(theta),
     ylab = "Kurtosis",
     main = "Kurtosis",
     ylim = c(0, 50))   

grid()
par(mfrow = c(1,1))



