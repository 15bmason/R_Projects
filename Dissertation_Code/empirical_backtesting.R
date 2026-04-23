library(rugarch)
library(quantmod)
library(moments)
library(tseries)

symbol <- "BTC-USD" 
getSymbols(symbol, src = "yahoo", from = "2018-01-01", to = "2026-01-01")

# Extract Adjusted Close and Calculate Log Returns
prices <- Cl(get(symbol))
real_returns <- diff(log(prices))[-1] # Removes first NA
real_returns <- as.numeric(real_returns)
total_obs <- length(real_returns)

plot(real_returns, type='l', main=paste("Log Returns of", symbol))
# Squared returns showing clustering
plot(real_returns^2, type='l', main="Volatility Clustering (Squared Returns)")

par(mfrow = c(2, 1), mar = c(4, 4, 2, 1)) 

# Log Returns (Panel A)
plot(real_returns, 
     type = "l", 
     main = "Bitcoin Daily Log Returns", 
     ylab = "Log Return", 
     xlab = "",
     col = "black", 
     lwd = 0.8)

# Squared Returns (Panel B)
plot(real_returns^2, 
     type = "l", 
     main = "Bitcoin Squared Daily Returns", 
     ylab = "Squared Return", 
     xlab = "Time (Days)", 
     col = "black",
     lwd = 0.8)

par(mfrow = c(1, 1))

# Forecasting (rolling window)
window_size <- 922
forecasts_linear <- numeric(total_obs)
forecasts_qmle   <- numeric(total_obs)

# Specification for rugarch
spec <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,0)), 
                   mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
                   distribution.model="norm")

cat("\nStarting Rolling Forecast\n")
for (t in window_size:(total_obs - 1)) {
  
  # Define the Data Window
  window_data <- real_returns[(t - window_size + 1):t]
  
  # Linear estimator forecast
  est <- linear_arch_estimator(window_data)
  beta0_lin <- est$WLS[1]
  beta1_lin <- est$WLS[2]
  
  current_return_sq <- real_returns[t]^2
  
  if(!is.na(beta0_lin)) {
    var_pred <- beta0_lin + beta1_lin * current_return_sq
    forecasts_linear[t+1] <- ifelse(var_pred > 0, var_pred, 0.00001)
  } else {
    forecasts_linear[t+1] <- NA
  }
  
  # QMLE Forecast
  try({
    fit <- ugarchfit(spec=spec, data=window_data, solver="hybrid")
    forc <- ugarchforecast(fit, n.ahead=1)
    forecasts_qmle[t+1] <- sigma(forc)^2 # Extract variance
  }, silent=TRUE)
  
  if(t %% 50 == 0) cat(t, "/", total_obs, "\n")
}

# Align data
test_data <- real_returns[(window_size + 1):total_obs]
sigma_linear <- sqrt(forecasts_linear[(window_size + 1):total_obs])
sigma_qmle   <- sqrt(forecasts_qmle[(window_size + 1):total_obs])

alpha <- 0.01
z_score <- qnorm(alpha)

VaR_linear <- z_score * sigma_linear
VaR_qmle   <- z_score * sigma_qmle

# Count violations
violations_linear <- sum(test_data < VaR_linear, na.rm=TRUE)
violations_qmle   <- sum(test_data < VaR_qmle, na.rm=TRUE)

cat("\n--- VaR Backtest Results (99% Confidence) ---\n")
cat("Total Days Forecasted:", length(test_data), "\n")
cat("Target Violations (1%):", round(length(test_data)*0.01), "\n\n")

# Count and percentage
cat("Linear Estimator Violations:", violations_linear, 
    "(", round(violations_linear/length(test_data)*100, 2), "%)\n")

cat("QMLE Estimator Violations:  ", violations_qmle, 
    "(", round(violations_qmle/length(test_data)*100, 2), "%)\n")

# Setup
y_limits <- range(c(test_data, VaR_linear, VaR_qmle), na.rm = TRUE)

# Returns
plot(test_data, type='l', col="darkgrey", 
     main="VaR: Linear (Red) and QMLE (Blue)", 
     ylab="Log Returns", xlab="Time", ylim=y_limits)

# QMLE VaR
lines(VaR_qmle, col="blue", lwd=1.5)

# Linear estimator VaR
lines(VaR_linear, col="red", lwd=1.5)

legend("bottomleft", 
       legend=c("Actual Returns", "Linear Est VaR", "QMLE VaR"),
       col=c("darkgrey", "red", "blue"), 
       lty=1, lwd=c(1, 1.5, 1.5), 
       bg="white", cex=0.8)

par(mfrow = c(2, 1), mar = c(2, 4, 2, 1), oma = c(2, 0, 2, 0))

# QMLE
plot(test_data, type = "l", col = "darkgrey", 
     main = "1% VaR Forecasts: QMLE vs. Linear WLS", 
     ylab = "Log Returns", xlab = "", 
     ylim = c(min(test_data) - 0.05, max(test_data) + 0.05))

lines(VaR_qmle, col = "blue", lwd = 1.5)

legend("topleft", legend = c("Actual Returns", "QMLE VaR"), 
       col = c("darkgrey", "blue"), lty = 1, lwd = 1.5, 
       bty = "n") 

# Linear WLS
plot(test_data, type = "l", col = "darkgrey", 
     main = "", 
     ylab = "Log Returns", xlab = "Time (Days)", 
     ylim = c(min(test_data) - 0.05, max(test_data) + 0.05))

lines(VaR_linear, col = "red", lwd = 1.5)

legend("topleft", legend = c("Actual Returns", "Linear WLS VaR"), 
       col = c("darkgrey", "red"), lty = 1, lwd = 1.5, 
       bty = "n")

par(mfrow = c(1, 1))

# Standard QMLE GARCH(1,1) Model
garch11_spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "norm"
)

# Fit the model to the full Bitcoin dataset
full_fit_garch11 <- ugarchfit(spec = garch11_spec, data = real_returns)
show(full_fit_garch11)

# Run the Rolling Window Backtest
roll_garch11 <- ugarchroll(
  spec = garch11_spec, 
  data = real_returns, 
  n.ahead = 1, 
  forecast.length = 2000,
  refit.every = 1, 
  refit.window = "moving", 
  window.size = 922,
  calculate.VaR = TRUE,
  VaR.alpha = 0.01   # 1% VaR
)

# Extract the VaR violation rate
report(roll_garch11, type = "VaR", VaR.alpha = 0.01, conf.level = 0.95)

# Plot used for presentation

total_days <- length(test_data)
start_plot <- max(1, total_days - 300) 
plot_idx   <- start_plot:total_days

data_recent <- test_data[plot_idx]
lin_recent  <- VaR_linear[plot_idx]
qmle_recent <- VaR_qmle[plot_idx]

par(mar = c(3, 4, 3, 1))

y_lims <- range(c(data_recent, lin_recent, qmle_recent), na.rm = TRUE)

plot(data_recent, type = 'l', 
     col = "darkgrey", lwd = 1.5, 
     ylim = y_lims,
     main = "Recent Performance (Last 300 Days)",
     ylab = "Log Returns", xlab = "Time",
     cex.main = 1.5, cex.lab = 1.2, cex.axis = 1.2)

# Add QMLE (Blue)
lines(qmle_recent, col = "blue", lwd = 1.5)

# Add Linear (Red)
lines(lin_recent, col = "red", lwd = 1.5)

# Add Legend
legend("bottomright", 
       legend = c("Market Returns", "Linear Est", "QMLE"),
       col = c("darkgrey", "red", "blue"), 
       lty = c(1, 1, 1), lwd = c(1.5, 1.5, 1.5), 
       bg = "white", cex = 1.1)

kupiec <- function(x, n, p) {
  p_hat <- x / n
  
  log_L_null <- (n - x) * log(1 - p) + x * log(p)
  
  if (x == 0) {
    log_L_alt <- n * log(1)
  } else if (x == n) {
    log_L_alt <- n * log(1)
  } else {
    log_L_alt <- (n - x) * log(1 - p_hat) + x * log(p_hat)
  }
  
  # Likelihood Ratio (LR) Statistic
  lr_stat <- 2 * (log_L_alt - log_L_null)
  
  # P-value from Chi-square distribution with 1 degree of freedom
  p_val <- 1 - pchisq(lr_stat, df = 1)
  
  return(list(stat = lr_stat, p = p_val))
}

alpha <- 0.01
kup_qmle <- kupiec(violations_qmle, length(test_data), alpha)
kup_le <- kupiec(violations_linear, length(test_data), alpha)

cat("QMLE LR Statistic: ", round(kup_qmle$stat, 4), "\n")
cat("QMLE P-Value: ", round(kup_qmle$p, 5), "\n")

cat("WLS LR Statistic: ", round(kup_le$stat, 4), "\n")
cat("WLS P-Value: ", round(kup_le$p, 5), "\n")
