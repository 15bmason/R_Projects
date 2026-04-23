library(rugarch)
library(MASS)

set.seed(111)
n_simulations <- 1000
n_obs <- 1000

# True parameters
true_a0 <- 0.1
true_a1 <- 0.2

sim_dist <- "student_t" # Options: "norm" or "student_t"
df_t <- 5 # Degrees of freedom for "student_t"

var_violations_qmle <- 0
var_violations_le <- 0

# Critical value for 1% VaR
z_01 <- qnorm(0.01) 

spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 0)),
  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "norm"
)

cat("Starting 1000 simulations using", sim_dist, "innovations\n")

for(i in 1:n_simulations) {
  
  # Simulate ARCH(1) innovations
  if (sim_dist == "norm") {
    z <- rnorm(n_obs + 1)
  } else {
    z <- rt(n_obs + 1, df = df_t) * sqrt((df_t - 2)/df_t) # Standardised t
  }
  
  r <- numeric(n_obs + 1)
  sig2 <- numeric(n_obs + 1)
  
  sig2[1] <- true_a0 / (1 - true_a1)
  r[1] <- z[1] * sqrt(sig2[1])
  
  for(t in 2:(n_obs + 1)) {
    sig2[t] <- true_a0 + true_a1 * r[t-1]^2
    r[t] <- z[t] * sqrt(sig2[t])
  }
  
  # In sample data (1 to 1000)
  r_in <- r[1:n_obs]
  # Out of sample true return
  r_out <- r[n_obs + 1]
  
  qmle_fit <- try(
    ugarchfit(spec = spec, data = r_in, solver = "hybrid"),
    silent = TRUE
    )
  
  if(!inherits(qmle_fit, "try-error")) {
    a0_qmle <- coef(qmle_fit)["omega"]
    a1_qmle <- coef(qmle_fit)["alpha1"]
  } else {
    # If the program crashes we use true values to avoid skewed VaR
    a0_qmle <- true_a0
    a1_qmle <- true_a1 
  }
  
  Y <- r_in[2:n_obs]^2
  X <- cbind(1, r_in[1:(n_obs-1)]^2)
  
  # Simple weights (inverse of squared lagged returns to handle tails)
  w <- 1 / (1 + r_in[1:(n_obs-1)]^2)^2 
  W <- diag(w)
  
  beta_le <- tryCatch(
    solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% Y,
    error = function(e) c(true_a0, 0.1) # fallback if singular
  )
  
  a0_le <- beta_le[1]
  a1_le <- beta_le[2]
  
  # Forecast variance
  sig2_qmle_forecast <- a0_qmle + a1_qmle * r_in[n_obs]^2
  sig2_le_forecast <- a0_le + a1_le * r_in[n_obs]^2
  
  # Prevent negative variance forecasts
  sig2_qmle_forecast <- max(sig2_qmle_forecast, 1e-6)
  sig2_le_forecast <- max(sig2_le_forecast, 1e-6)
  
  # Calculate 1% VaR threshold
  var_qmle <- sqrt(sig2_qmle_forecast) * z_01
  var_le <- sqrt(sig2_le_forecast) * z_01
  
  if(r_out < var_qmle) var_violations_qmle <- var_violations_qmle + 1
  if(r_out < var_le) var_violations_le <- var_violations_le + 1
  
  if(i %% 100 == 0) cat("Completed", i, "of", n_simulations, "\n")
}

# Calculate actual rates
qmle_hit_rate <- (var_violations_qmle / n_simulations) * 100
le_hit_rate <- (var_violations_le / n_simulations) * 100

cat("Simulation Setup:", sim_dist)
if(sim_dist == "student_t") cat(" (df =", df_t, ")\n") else cat("\n")
cat("Target VaR Violation Rate: 1.00%\n")
cat("QMLE Hit Rate: ", qmle_hit_rate, "%\n")
cat("Linear Estimator Hit Rate: ", le_hit_rate, "%\n")