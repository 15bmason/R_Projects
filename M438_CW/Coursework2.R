sim_birth_death_process <- function(n_0, fn_birth_rate, params_birth,
                                    fn_death_rate, params_death, end_time){
  n <- n_0
  time <- 0
  i <- 1
  t_vec <- numeric(0)
  X <- integer(0)
  t_vec[i] <- time
  X[i] <- n
  
  while (n < end_time){
    # Calculate birth and death rates for current value of n
    lambda_n = fn_birth_rate(n, params_birth)
    mu_n = fn_death_rate(n, params_death)
    # Simulate time until next event
    t_next <- rexp(1, lambda_n + mu_n)
    # Determine whether next event is a birth or a death
    u <- runif(1)
    if (lambda_n/(lambda_n + mu_n) > u){ # birth
      n <- n + 1
    } else { # death
      n <- n - 1
    }
    # Update time
    time <- time + t_next
    # Store time of event
    t_vec[i + 1] <- time
    # Store value of process
    X[i + 1] <- n
    # Increment counter
    i <- i + 1
  }
  return(list(X = matrix(X[-length(X)], nrow = 1), t = t_vec[-length(t_vec)], final_time = time)) 
}

# Set the birth rate and death rate functions to be those for the linear
# birth-death process
birth_rate <- function(n, params_birth){
  return(params_birth$lambda * (n - sqrt(n)))
}

death_rate <- function(n, params_death){
  return(params_death$mu * n)
}

n_0 <- 5
end_time <- 1000
params_birth <- list(lambda = 1)
params_death <- list(mu = 0)
set.seed(2)

n_sims <- 10000
sim_T <- numeric(n_sims)

for (j in 1:n_sims) {
  res_bd <- sim_birth_death_process(n_0, birth_rate, params_birth,
                                    death_rate, params_death, end_time)
  sim_T[j] <- res_bd$final_time
}

cat("Mean:", mean(sim_T), "\n")
cat("Var:", var(sim_T), "\n")

# part c
cat(sum(1 / (1 * (5:999 - sqrt(5:999)))))
cat(sum(1 / (1 * (5:999 - sqrt(5:999))^2)))

