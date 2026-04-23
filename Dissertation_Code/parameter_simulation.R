library(rugarch)
library(quantmod)
library(moments)
library(tseries)

linear_arch_estimator <- function(x_series) {
  n <- length(x_series)
  sq_returns <- x_series^2
  Y_vec <- sq_returns[2:n]
  
  # Col 1 is intercept, Col 2 is lagged squared return
  Z_mat <- cbind(1, sq_returns[1:(n-1)])
  
  # OLS
  bhat_pr <- tryCatch({
    solve(crossprod(Z_mat), crossprod(Z_mat, Y_vec))
  }, error = function(e) return(c(NA, NA)))
  
  if (any(is.na(bhat_pr))) return(list(OLS = c(NA, NA), WLS = c(NA, NA)))
  
  # FGLS
  estimated_var <- Z_mat %*% bhat_pr
  
  estimated_var[estimated_var <= 0] <- 1e-6 # Error prevention
  
  Z_weighted <- Z_mat / as.vector(estimated_var)
  Y_weighted <- Y_vec / as.vector(estimated_var)
  
  bhat_wls <- tryCatch({
    solve(crossprod(Z_weighted), crossprod(Z_weighted, Y_weighted))
  }, error = function(e) return(c(NA, NA)))
  
  return(list(OLS = as.vector(bhat_pr), WLS = as.vector(bhat_wls)))
}

# Simulation parameters
n_sims <- 1000
n_obs  <- 1000
beta0_true <- 0.1
beta1_true <- 0.2

# Storage
results_linear_ols <- matrix(NA, nrow = n_sims, ncol = 2)
results_linear_wls <- matrix(NA, nrow = n_sims, ncol = 2)
results_qmle       <- matrix(NA, nrow = n_sims, ncol = 2)

colnames(results_linear_ols) <- c("beta0", "beta1")
colnames(results_linear_wls) <- c("beta0", "beta1")
colnames(results_qmle)       <- c("beta0", "beta1")

# QMLE
spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 0)),
  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "norm"
)

set.seed(123)
for(i in 1:n_sims) {
  
  # DISTRIBUTIONS --------------------------------------------------------------
  
  # Normal =====================================================================
  #epsilon <- rnorm(n_obs)  # Normal
  
  # Student's t/Heavy tails ====================================================
  #epsilon <- rt(n_obs, df=5) * sqrt(3/5)
  
  # Logistic ===================================================================
  #epsilon <- rlogis(n_obs) * (sqrt(3)/pi)
  
  # Double exponential =========================================================
  u <- runif(n_obs) - 0.5 
  epsilon <- -sign(u) * log(1 - 2*abs(u)) * (1/sqrt(2))
  
  # ----------------------------------------------------------------------------

  X <- numeric(n_obs)
  sigma <- numeric(n_obs)
  
  sigma[1] <- sqrt(beta0_true / (1 - beta1_true))
  X[1] <- sigma[1] * epsilon[1]
  
  for(t in 2:n_obs){
    sigma[t] <- sqrt(beta0_true + beta1_true * (X[t-1])^2)
    X[t] <- sigma[t] * epsilon[t]
  }
  
  # Linear estimation
  est <- linear_arch_estimator(X)
  results_linear_ols[i, ] <- est$OLS
  results_linear_wls[i, ] <- est$WLS
  
  # QMLE - try preventing any crashes
  try({
    fit <- ugarchfit(spec = spec, data = X, solver = "hybrid")
    coefs <- coef(fit)
    results_qmle[i, ] <- c(coefs["omega"], coefs["alpha1"])
  }, silent = TRUE)
  
}

# Computing metrics function
calc_metrics <- function(estimates, true_val) {
  clean_est <- na.omit(estimates)
  
  mean_est <- mean(clean_est)
  bias     <- mean_est - true_val
  variance <- var(clean_est)
  mse      <- variance + bias^2
  
  return(c(Mean = mean_est, Bias = bias, MSE = mse))
}

cat("\n\n--- SIMULATION RESULTS ---\n")
cat("Sample Size:", n_obs, "| Sims:", n_sims, "\n")
cat("True Beta1 (Alpha1):", beta1_true, "\n\n")

b1_ols  <- calc_metrics(results_linear_ols[, 2], beta1_true)
b1_wls  <- calc_metrics(results_linear_wls[, 2], beta1_true)
b1_qmle <- calc_metrics(results_qmle[, 2], beta1_true)

# We look at Column 1 of the results matrices
b0_ols  <- calc_metrics(results_linear_ols[, 1], beta0_true)
b0_wls  <- calc_metrics(results_linear_wls[, 1], beta0_true)
b0_qmle <- calc_metrics(results_qmle[, 1], beta0_true)

# Create the Table for Beta 0
table_b0 <- rbind(OLS = b0_ols, WLS = b0_wls, QMLE = b0_qmle)
table_b1 <- rbind(OLS = b1_ols, WLS = b1_wls, QMLE = b1_qmle)

cat("\n RESULTS FOR BETA 0 (Intercept) \n")
print(table_b0)

cat("\n RESULTS FOR BETA 1 (ARCH Parameter) \n")
print(table_b1)

# BOXPLOTS 
par(mfrow=c(1,2))

# Beta 0
boxplot(results_linear_ols[,1], results_linear_wls[,1], results_qmle[,1],
        names = c("OLS", "WLS", "QMLE"),
        main = expression(paste("Intercept ", alpha[0])),
        col = c("lightblue", "lightgreen", "salmon"))
abline(h = beta0_true, col = "red", lty = 2)

# Beta 1
boxplot(results_linear_ols[,2], results_linear_wls[,2], results_qmle[,2],
        names = c("OLS", "WLS", "QMLE"),
        main = expression(paste("ARCH Parameter ", alpha[1])),
        col = c("lightblue", "lightgreen", "salmon"))
abline(h = beta1_true, col = "red", lty = 2)

par(mfrow=c(1,1))

boxplot(results_linear_ols[,2], results_linear_wls[,2], results_qmle[,2],
        names = c("OLS", "WLS", "QMLE"),
        main = expression(paste("ARCH Parameter ", alpha[1])),
        ylab = "Estimated Value",
        col = c("lightblue", "lightgreen", "salmon"))
abline(h = beta1_true, col = "red", lty = 2, lwd = 2) # True value line

linear_arch_estimator <- function(x_series) {
  n <- length(x_series)
  
  sq_returns <- x_series^2
  
  
  Y_vec <- sq_returns[2:n]
  
  
  Z_mat <- cbind(1, sq_returns[1:(n-1)])
  
  # OLS
  bhat_pr <- tryCatch({
    solve(crossprod(Z_mat), crossprod(Z_mat, Y_vec))
  }, error = function(e) return(c(NA, NA)))
  
  if (any(is.na(bhat_pr))) {
    return(list(OLS = c(NA, NA), WLS = c(NA, NA)))
  }
  
  # WLS
  estimated_var <- Z_mat %*% bhat_pr
  estimated_var[estimated_var <= 0] <- 1e-6
  
  Z_weighted <- Z_mat / as.vector(estimated_var)
  Y_weighted <- Y_vec / as.vector(estimated_var)
  
  bhat_wls <- tryCatch({
    solve(crossprod(Z_weighted), crossprod(Z_weighted, Y_weighted))
  }, error = function(e) {
    return(c(NA, NA))
  })
  
  return(list(OLS = as.vector(bhat_pr), WLS = as.vector(bhat_wls)))
}