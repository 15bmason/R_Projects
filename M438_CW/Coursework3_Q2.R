# Generative AI was used as a general helping tool for this question

fn_rate <- function(X, time, params){
  P2 <- X[1]
  c_val <- params$c
  
  rate_fwd <- params$k1 * max(0, c_val - 2*P2) * max(0, c_val - 2*P2 - 1) / 2
  rate_rev <- params$k2 * P2
  
  return(c(rate_fwd, rate_rev))
}

S <- matrix(c(1, -1), nrow = 1, ncol = 2)
fn_update_state <- make_update_state(S)

params <- list(k1 = 1.66e-3, k2 = 0.2, c = 301)
X_0 <- matrix(0, nrow = 1, ncol = 1)
end_time <- 15
n_sims <- 1000

final_P2 <- numeric(n_sims)
trajectories <- vector("list", n_sims)

set.seed(123456)
cat("Running", n_sims, "simulations.\n")

for(i in 1:n_sims){
  res <- SSA(X_0, fn_rate, fn_update_state, make_sim_one_event, params, end_time)
  trajectories[[i]] <- res
  final_P2[i] <- res$X[length(res$X)] 
}

est_prob <- mean(final_P2 >= 70 & final_P2 <= 90)
est_mean <- mean(final_P2)

cat(sprintf("Estimated P(70 <= P2(15) <= 90): %.3f\n", est_prob))
cat(sprintf("Estimated E[P2(15)]: %.3f\n", est_mean))

max_P2 <- floor(params$c / 2) # 150
n_vals <- 0:max_P2

log_w <- numeric(length(n_vals))
log_w[1] <- 0 

for (i in 1:max_P2) {
  n <- i - 1
  log_w[i + 1] <- log_w[i] + log(params$k1) + log(params$c - 2*n) + 
    log(params$c - 2*n - 1) - log(2) - log(params$k2) - log(i)
}

log_w <- log_w - max(log_w)
w <- exp(log_w)
pmf <- w / sum(w)

theo_mean <- sum(n_vals * pmf)
theo_var <- sum((n_vals - theo_mean)^2 * pmf)
theo_sd <- sqrt(theo_var)

plot(0, 0, type = "n", xlim = c(0, end_time), ylim = c(0, 110), 
     xlab = "Time (s)", ylab = expression(P[2](t)), 
     main = "Simulated Trajectories")

for(i in 1:n_sims){
  lines(trajectories[[i]]$t, trajectories[[i]]$X, type = "s", col = rgb(0, 0, 0, 0.05))
}

abline(h = theo_mean, col = "red", lwd = 2)
abline(h = theo_mean + 3*theo_sd, col = "red", lty = 2, lwd = 2)
abline(h = theo_mean - 3*theo_sd, col = "red", lty = 2, lwd = 2)

legend("bottomright", legend = c("Trajectories", "Theo Mean", "+/- 3 SD"), 
       col = c("gray", "red", "red"), lty = c(1, 1, 2), lwd = c(1, 2, 2), cex=0.7)

hist(final_P2, breaks = seq(-0.5, max_P2 + 0.5, by = 1), probability = TRUE,
     xlim = c(50, 110), col = "lightblue", border = "white",
     xlab = expression(P[2](15)), main = "Distribution at T=15")

lines(n_vals, pmf, type = "o", col = "red", pch = 16, cex = 0.5)

legend("topright", legend = c("Simulated", "Theoretical PMF"),
       fill = c("lightblue", NA), border = c("white", NA),
       lty = c(NA, 1), pch = c(NA, 16), col = c(NA, "red"), cex = 0.7)
