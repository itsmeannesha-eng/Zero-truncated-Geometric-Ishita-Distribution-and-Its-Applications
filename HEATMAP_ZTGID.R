library(ggplot2)
library(reshape2)

# Mean function
mu_zgid <- function(theta){
  (theta^3/(theta^3 + 2)) * ( theta/(theta - 1) + 2/(theta - 1)^3 )
}

# Variance function
var_ztgid <- function(theta){
  
  term1 <- (2*theta^3/(theta^3+2)) * ( 2/(theta-2)^3 + theta/(theta-2) )
  
  term2 <- (theta^3/(theta^3+2)) * ( 2/(theta-1)^3 + theta/(theta-1) )
  
  term3 <- (theta^6/(theta^3+2)^2) * ( theta/(theta-1) + 2/(theta-1)^3 )^2
  
  term1 - term2 - term3
}

# theta values
theta <- seq(2.1,6,length=200)

# compute mean and variance
mu <- mu_zgid(theta)
variance <- var_ztgid(theta)

# create dataframe
df <- data.frame(
  theta = theta,
  Mean = mu,
  Variance = variance
)

# reshape
df_long <- melt(df, id.vars="theta")

# extreme values
cap <- quantile(df_long$value, 0.95, na.rm=TRUE)
df_long$value_cap <- pmin(df_long$value, cap)

# heatmap
ggplot(df_long, aes(x = theta, y = variable, fill = value_cap)) +
  geom_tile() +
  scale_fill_gradient(low = "red", high = "green") +
  labs(x = expression(theta),
       y = "",
       fill = "Value",
       title = "Heatmap of Mean and Variance of ZTGID") +
  theme_minimal()

