library(quantmod)
library(rugarch)

# Fetch the data
getSymbols("BTC-USD", src = "yahoo", from = "2018-01-01", to = "2026-01-01")
prices <- Cl(get("BTC-USD"))
real_returns <- as.numeric(diff(log(prices))[-1])
n <- length(real_returns)

mean_returns <- mean(real_returns)
sd_returns <- sd(real_returns)
min_returns <- min(real_returns)
max_returns <- max(real_returns)
skew_returns <- skewness(real_returns)
kurt_returns <- kurtosis(real_returns)

cat("Observations (N):", n, "\n")
cat("Mean:", mean_returns, "\n")
cat("Std Dev:", sd_returns, "\n")
cat("Minimum:", min_returns, "\n")
cat("Maximum:", max_returns, "\n")
cat("Skewness:", skew_returns, "\n")
cat("Kurtosis:", kurt_returns, "\n")

# Null Hypothesis: Data is non stationary
cat("\n--- ADF Test ---\n")
print(adf.test(real_returns))

# Null Hypothesis: No autocorrelation in variance
cat("\n--- Ljung-Box Test (Squared Returns) ---\n")
print(Box.test(real_returns^2, lag = 20, type = "Ljung-Box"))

cat("Total observations:", n, "\n\n")

spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 0)), 
                   mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                   distribution.model = "norm")

# Fit QMLE on the entire dataset
fit_qmle <- ugarchfit(spec = spec, data = real_returns, solver = "hybrid")

# Extract QMLE Coefficients and Standard Errors
qmle_matcoef <- fit_qmle@fit$robust.matcoef
qmle_alpha0_est <- qmle_matcoef["omega", 1]
qmle_alpha0_se  <- qmle_matcoef["omega", 2]
qmle_alpha1_est <- qmle_matcoef["alpha1", 1]
qmle_alpha1_se  <- qmle_matcoef["alpha1", 2]

sq_returns <- real_returns^2
Y_vec <- sq_returns[2:n]
Z_mat <- cbind(1, sq_returns[1:(n-1)])

# OLS
bhat_ols <- solve(crossprod(Z_mat), crossprod(Z_mat, Y_vec))

# WLS
estimated_var <- Z_mat %*% bhat_ols
estimated_var[estimated_var <= 0] <- 1e-6 # Enforce positivity

# Divide Z and Y by variance creates the 1/sigma^4 weight
Z_weighted <- Z_mat / as.vector(estimated_var)
Y_weighted <- Y_vec / as.vector(estimated_var)

# Final WLS Estimates
bhat_wls <- solve(crossprod(Z_weighted), crossprod(Z_weighted, Y_weighted))

# Calculate the residuals
u_hat <- Y_vec - (Z_mat %*% bhat_wls)

# Standardise residuals
sigma2_hat <- Z_mat %*% bhat_wls
sigma2_hat[sigma2_hat <= 0] <- 1e-6
eta_hat <- u_hat / sigma2_hat

# Variance of eta
var_eta <- var(eta_hat)

# Covariance Matrix
cov_matrix <- as.numeric(var_eta) * solve(crossprod(Z_weighted))

# Standard Errors are the square root of the diagonal
wls_se <- sqrt(diag(cov_matrix))

wls_alpha0_est <- bhat_wls[1]
wls_alpha0_se  <- wls_se[1]
wls_alpha1_est <- bhat_wls[2]
wls_alpha1_se  <- wls_se[2]

cat(" FINAL PARAMETER ESTIMATES \n")
cat("QMLE (Standard)\n")
cat(sprintf("  Alpha 0 (Intercept): %.6f  (SE: %.6f)\n", qmle_alpha0_est, qmle_alpha0_se))
cat(sprintf("  Alpha 1 (ARCH Term): %.6f  (SE: %.6f)\n\n", qmle_alpha1_est, qmle_alpha1_se))

cat("Linear WLS (Robust)\n")
cat(sprintf("  Alpha 0 (Intercept): %.6f  (SE: %.6f)\n", wls_alpha0_est, wls_alpha0_se))
cat(sprintf("  Alpha 1 (ARCH Term): %.6f  (SE: %.6f)\n", wls_alpha1_est, wls_alpha1_se))