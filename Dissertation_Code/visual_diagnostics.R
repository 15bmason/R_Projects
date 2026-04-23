library(quantmod)

# Fetch data
getSymbols("BTC-USD", src = "yahoo", from = "2018-01-01", to = "2026-01-01")
prices <- Cl(get("BTC-USD"))
real_returns <- as.numeric(diff(log(prices))[-1])

cat("Q-Q Plot \n")

# Standardise the returns to compare against a N(0,1) distribution
std_returns <- (real_returns - mean(real_returns)) / sd(real_returns)

qqnorm(std_returns, 
       main = "Normal Q-Q Plot: Bitcoin Log Returns",
       xlab = "Theoretical Normal Quantiles",
       ylab = "Empirical Bitcoin Quantiles",
       pch = 20, 
       col = "blue") 
qqline(std_returns, col = "red", lwd = 2)
grid()
