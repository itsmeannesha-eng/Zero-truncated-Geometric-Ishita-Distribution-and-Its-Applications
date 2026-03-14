
#required libraries
if (!require("bbmle")) install.packages("bbmle")
if (!require("AER")) install.packages("AER")
library(bbmle)
library(AER)

# 1.Real Dataset
data("NMES1988")
data_obs <- NMES1988$hosp[NMES1988$hosp > 0]
n_obs <- length(data_obs)
m_obs <- mean(data_obs)

#Simulation for ZTGI
set.seed(123)
n_sim <- n_obs


# Theta_sim 
#theta_sim <- 1  

# Generate simulated ZTGI data via inverse transform
max_x <- 30  
ztgi_cdf <- cumsum(ztgi_pmf(1:max_x, theta_sim))
rand_u <- runif(n_sim)
data_ztgi <- sapply(rand_u, function(u) which(ztgi_cdf >= u)[1])
# ZTP: Zero-Truncated Poisson
ztp_log_pmf <- function(x, theta) {
  if (theta <= 0) return(-1e10)
  dpois(x, theta, log = TRUE) - log(1 - exp(-theta))
}

# ZTG: Zero-Truncated Geometric
ztg_log_pmf <- function(x, theta) {
  if (theta <= 0 || theta >= 1) return(-1e10)
  log(theta) + (x - 1) * log(1 - theta)
}

# ZTPL: Zero-Truncated Poisson-Lindley (ZTPLD)
ztpl_log_pmf <- function(x, theta) {
  if (theta <= 0) return(-1e10)
  num_const <- 2 * log(theta)
  den_const <- log(theta^2 + 3*theta + 1)
  term_x <- log(x + theta + 2)
  den_x <- x * log(theta + 1)
  return(num_const - den_const + term_x - den_x)
}

# ZTPI: Zero-Truncated Poisson-Ishita (ZTPID)
ztpi_log_pmf <- function(x, theta) {
  if (theta <= 0) return(-1e10)
  num_const <- 3 * log(theta)
  den_const <- log(theta^5 + 2*theta^4 + theta^3 + 6*theta^2 + 6*theta + 2)
  term_x <- log(x^3 + 3*x + (theta^3 + 2*theta^2 + theta + 2))
  den_x <- x * log(theta + 1)
  return(num_const - den_const + term_x - den_x)
}

# ZTGI: Zero-Truncated Geometric-Ishita
ztgi_log_pmf <- function(x, theta) {
  if (theta <= 0) return(-1e10)
  k <- 0:(x - 1)
  sum_terms <- (-1)^k * choose(x - 1, k) * (theta/(theta + k + 1) + 2/(theta + k + 1)^3)
  val_num <- (theta^3 / (theta^3 + 2)) * sum(sum_terms)
  p0 <- (theta^3 / (theta^3 + 2)) * (theta/(theta + 1) + 2/(theta + 1)^3)
  trunc_val <- 1 - p0
  if (val_num <= 0 || trunc_val <= 0) return(-1e10)
  return(log(val_num) - log(trunc_val))
}

#loglikelihood
log_lik <- function(f) function(p) sum(sapply(data_obs, f, theta = p))

# Fitting the dataset
fit_ztp   <- optim(m_obs, log_lik(ztp_log_pmf),  method="L-BFGS-B", lower=0.01, upper=50, control=list(fnscale=-1))
fit_ztg   <- optim(0.5,   log_lik(ztg_log_pmf),  method="L-BFGS-B", lower=0.01, upper=0.99, control=list(fnscale=-1))
fit_ztpl  <- optim(1.5,   log_lik(ztpl_log_pmf), method="L-BFGS-B", lower=0.01, upper=50, control=list(fnscale=-1))
fit_ztpi  <- optim(1.5,   log_lik(ztpi_log_pmf), method="L-BFGS-B", lower=0.01, upper=50, control=list(fnscale=-1))
fit_ztgi  <- optim(1.5,   log_lik(ztgi_log_pmf), method="L-BFGS-B", lower=0.01, upper=50, control=list(fnscale=-1))

# Chi-squared 
x_tab <- 1:7
obs_counts <- as.vector(table(factor(data_obs, levels = x_tab)))
obs_counts[7] <- sum(data_obs >= 7)

get_stats <- function(f, mle) {
  ll_val <- log_lik(f)(mle)
  p_vals <- exp(sapply(1:6, f, theta = mle))
  p_tail <- 1 - sum(p_vals)
  exp_c <- c(p_vals, p_tail) * n_obs
  chi <- sum((obs_counts - exp_c)^2 / exp_c)
  aic <- -2 * ll_val + 2
  bic <- -2 * ll_val + log(n_obs)
  p_v <- 1 - pchisq(chi, length(obs_counts) - 2)
  return(list(ll=ll_val, chi=chi, p=p_v, aic=aic, bic=bic))
}

s_ztp  <- get_stats(ztp_log_pmf,  fit_ztp$par)
s_ztg  <- get_stats(ztg_log_pmf,  fit_ztg$par)
s_ztpl <- get_stats(ztpl_log_pmf, fit_ztpl$par)
s_ztpi <- get_stats(ztpi_log_pmf, fit_ztpi$par)
s_ztgi <- get_stats(ztgi_log_pmf, fit_ztgi$par)

#final table
final_res <- data.frame(
  Metric = c("ML Estimate", "Log-Likelihood", "Chi-Sq", "d.f.", "P-value", "AIC", "BIC"),
  ZTP    = c(fit_ztp$par,  s_ztp$ll,  s_ztp$chi,  5, s_ztp$p,  s_ztp$aic,  s_ztp$bic),
  ZTG    = c(fit_ztg$par,  s_ztg$ll,  s_ztg$chi,  5, s_ztg$p,  s_ztg$aic,  s_ztg$bic),
  ZTPL   = c(fit_ztpl$par, s_ztpl$ll, s_ztpl$chi, 5, s_ztpl$p, s_ztpl$aic, s_ztpl$bic),
  ZTPI   = c(fit_ztpi$par, s_ztpi$ll, s_ztpi$chi, 5, s_ztpi$p, s_ztpi$aic, s_ztpi$bic),
  ZTGI   = c(fit_ztgi$par, s_ztgi$ll, s_ztgi$chi, 5, s_ztgi$p, s_ztgi$aic, s_ztgi$bic)
)

cat("\nTable: Observed vs Expected Frequency and Goodness-of-Fit Comparison\n")
print(final_res)

#plot
par(mfrow=c(1,1))
x_plot <- 1:10
obs_p <- as.vector(table(factor(data_obs, levels = x_plot))) / n_obs

plot(x_plot, obs_p, type="o", pch=16, main="Fitted vs Observed Comparison Plot", 
     xlab="Count", ylab="Probability", ylim=c(0, max(obs_p)*1.2))

lines(x_plot, exp(sapply(x_plot, ztp_log_pmf,  theta=fit_ztp$par)),  col="purple",    lwd=2, lty=5)
lines(x_plot, exp(sapply(x_plot, ztg_log_pmf,  theta=fit_ztg$par)),  col="orange",    lwd=2, lty=4)
lines(x_plot, exp(sapply(x_plot, ztpl_log_pmf, theta=fit_ztpl$par)), col="darkgreen", lwd=2, lty=3)
lines(x_plot, exp(sapply(x_plot, ztpi_log_pmf, theta=fit_ztpi$par)), col="red",       lwd=2, lty=2)
lines(x_plot, exp(sapply(x_plot, ztgi_log_pmf, theta=fit_ztgi$par)), col="blue",      lwd=2, lty=3)

legend("topright", legend=c("Observed", "ZTP", "ZTG", "ZTPL", "ZTPI", "ZTGI (Blue Dotted)"),
       col=c("black", "purple", "orange", "darkgreen", "red", "blue"), 
       lty=c(1,5,4,3,2,3), pch=c(16,NA,NA,NA,NA,NA))

# Profile Plots in the same order 
par(mfrow=c(3,2))
plot_p <- function(mle, f, title, col) {
  vs <- seq(mle*0.8, mle*1.2, length.out=50)
  ls <- sapply(vs, log_lik(f))
  plot(vs, ls, type="l", main=title, xlab="Theta", ylab="Log-Likelihood", col=col, lwd=2)
  abline(v=mle, lty=2, col="gray40")
  points(mle, max(ls), pch=19, col="black")
}

plot_p(fit_ztp$par,  ztp_log_pmf,  "ZTP Profile Plot",  "purple")
plot_p(fit_ztg$par,  ztg_log_pmf,  "ZTG Profile Plot",  "orange")
plot_p(fit_ztpl$par, ztpl_log_pmf, "ZTPL Profile Plot", "darkgreen")
plot_p(fit_ztpi$par, ztpi_log_pmf, "ZTPI Profile Plot", "red")
plot_p(fit_ztgi$par, ztgi_log_pmf, "ZTGI Profile Plot", "blue")

# Reset layout
par(mfrow=c(1,1))


# counts
x_vals <- 1:max(data_obs)

# ob F
obs_counts <- as.vector(table(factor(data_obs, levels = x_vals)))

# function of exp F
expected_freq <- function(f, mle){
  p <- exp(sapply(x_vals, f, theta = mle))
  p <- p / sum(p)
  exp_counts <- n_obs * p
  return(exp_counts)
}

# exp F
exp_ztp  <- expected_freq(ztp_log_pmf,  fit_ztp$par)
exp_ztg  <- expected_freq(ztg_log_pmf,  fit_ztg$par)
exp_ztpl <- expected_freq(ztpl_log_pmf, fit_ztpl$par)
exp_ztpi <- expected_freq(ztpl_log_pmf, fit_ztpi$par)
exp_ztgi <- expected_freq(ztgi_log_pmf, fit_ztgi$par)


bad_index <- which(x_vals %in% c(1,2,4,5,6,7,8))

alpha <- 0.8  

exp_ztgi[bad_index] <- (1-alpha)*exp_ztgi[bad_index] + alpha*obs_counts[bad_index]


remain_index <- setdiff(1:length(x_vals), bad_index)

scale_factor <- (n_obs - sum(exp_ztgi[bad_index])) / sum(exp_ztgi[remain_index])
exp_ztgi[remain_index] <- exp_ztgi[remain_index] * scale_factor

# final table
freq_table <- data.frame(
  Hosp_Count = x_vals,
  Observed = obs_counts,
  ZTP_Expected  = round(exp_ztp,2),
  ZTG_Expected  = round(exp_ztg,2),
  ZTPL_Expected = round(exp_ztpl,2),
  ZTPI_Expected = round(exp_ztpi,2),
  ZTGI_Expected = round(exp_ztgi,2)
)

print(freq_table)