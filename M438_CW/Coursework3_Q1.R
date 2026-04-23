# Use of generative AI to help with data handling

make_update_state <- function(S){
  fn_update_state <- function(X, time, j){
    X <- X + S[,j]
    return(X)
  }
  return(fn_update_state)
}

make_sim_one_event <- function(fn_rate, fn_update_state){
  fn_sim_one_event <- function(X, time, params){
    rate <- fn_rate(X, time, params)
    
    # Simulate time until next event
    u1 <- runif(1)
    t_next <- -log(u1)/sum(rate)
    
    # Determine what type of event next event is
    u2 <- runif(1)
    j <- min(which(u2 < cumsum(rate/sum(rate))))
    
    # Update state
    X <- fn_update_state(X, time, j)
    
    # Update time
    time <- time + t_next
    
    return(list(X = X, time = time, rate = rate))
  }
  return(fn_sim_one_event)
}

SSA <- function(X_0, fn_rate, fn_update_state, fn_make_sim_one_event, params, end_time){
  time <- 0
  i <- 1
  t_vec <- numeric(0)
  X_list <- vector("list", 0)
  
  t_vec[i] <- time
  X_list[[i]] <- X_0
  X <- X_0
  
  fn_sim_one_event <- fn_make_sim_one_event(fn_rate, fn_update_state)
  
  while (time < end_time && !all(fn_rate(X, time, params)==0)){
    res <- fn_sim_one_event(X, time, params)
    X <- res$X
    time <- res$time
    
    if (time < end_time){
      t_vec[i+1] <- time
      # Store value of state
      X_list[[i+1]] <- X
      i <- i+1
    }
  }
  X_arr <- do.call(cbind, X_list)
  return(list(X = X_arr, t = t_vec))
}

# PART A
fn_rate <- function(X, time, params){
  rate_birth <- params$lambda * X[1]
  rate_death <- params$mu * X[1]
  rate_immig <- params$nu
  return(c(rate_birth, rate_death, rate_immig))
}

S <- matrix(c(1, -1, 1), nrow = 1, ncol = 3)
fn_update_state <- make_update_state(S)

X_0 <- matrix(3, nrow = 1, ncol = 1) 
params <- list(lambda = 0.2, mu = 0.3, nu = 1)
end_time <- 10000

set.seed(123456)
sim_results <- SSA(X_0, fn_rate, fn_update_state, make_sim_one_event, params, end_time)

X_vec <- as.numeric(sim_results$X)

# PART B
max_x <- max(X_vec)
x_vals <- 0:max_x

r_size <- params$nu / params$lambda
p_prob <- 1 - (params$lambda / params$mu)
pi_theoretical <- dnbinom(x_vals, size = r_size, prob = p_prob)

hist(X_vec, breaks = seq(-0.5, max_x + 0.5, by = 1), probability = TRUE,
     main = "Part B: Initial Histogram", 
     xlab = "State (Population Size)", col = "lightblue", border = "white")
lines(x_vals, pi_theoretical, type = "o", col = "red", lwd = 2, pch = 16)
legend("topright", legend = c("Simulated", "Theoretical"),
       fill = c("lightblue", NA), border = c("white", NA),
       lty = c(NA, 1), pch = c(NA, 16), col = c(NA, "red"), cex = 0.8)

# PART C
t_vec <- sim_results$t
t_vec_extended <- c(t_vec, end_time)
waiting_times <- diff(t_vec_extended)

empirical_probs <- numeric(max_x + 1)
for (i in 0:max_x) {
  empirical_probs[i + 1] <- sum(waiting_times[X_vec == i]) / end_time
}

plot(x_vals, empirical_probs, type = "h", lwd = 5, col = "lightgreen",
     main = "Part C: Time Weighted", 
     xlab = "State (Population Size)", ylab = "Probability")
lines(x_vals, pi_theoretical, type = "o", col = "red", lwd = 2, pch = 16)
legend("topright", legend = c("Time Weighted", "Theoretical"),
       col = c("lightgreen", "red"), lwd = c(5, 2), pch = c(NA, 16), cex = 0.8)
