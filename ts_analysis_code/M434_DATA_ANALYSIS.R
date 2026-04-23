library(forecast)
library(tseries)
library(ggplot2)

# DATA LOADING =================================================================
raw_data <- read.csv("ei_eteu27_2020_m_page_linear_2_0.csv")
clean_data<-subset(raw_data,geo == "DE") # Filter data to find Germany
trade_values<-clean_data$OBS_VALUE # Pulls the trade values

full_ts_data <- ts(trade_values, start = c(2002, 1), frequency = 12)
ts_data <- window(full_ts_data, end = c(2024, 12)) # Cutoff at 2024
actual_2025 <- window(full_ts_data, start = c(2025, 1), end = c(2025, 3))
connector_line <- window(full_ts_data, start = c(2024, 12), end = c(2025, 1))
# Connector line is purely cosmetic so that there isn't a gap when plotting

# Plotting raw data
autoplot(full_ts_data) + 
  ggtitle("Trade Data 2002 - 2025") + 
  ylab("Trade Value (Million of Euro)") +
  theme_minimal(base_size = 16)

# BOX JENKINS ==================================================================
# Stationarity test

# Plot ACF and PACF Raw Series
ggtsdisplay(ts_data, main="Raw Data")
adf_raw <- adf.test(ts_data)
print(paste("ADF p-value:", round(adf_raw$p.value, 4)))

#autoplot(acf(ts_data)) +
  #theme_minimal(base_size=16)

# Differenced and Log Differenced Series
ts_diff <- diff(ts_data, differences = 1)          
ts_double_diff <- diff(ts_diff, lag = 12)

ts_log_diff <- diff(log(ts_data), differences=1)
ts_log_double_diff <- diff(ts_log_diff, lag=12)

# Plot Transformed and Differenced Series
autoplot(ts_log_double_diff) + 
  ggtitle("T&D Extra Import Trade Data 2002-2024") +
  ylab("T&D Trade Value (Millions of Euro)") + 
  theme_minimal(base_size=16)

# Plot ACF and PACF Single Differenced, ADF Test
ggtsdisplay(ts_diff, main="Single Differenced Data")
adf_single_diff <- adf.test(ts_diff)
print(paste("ADF p-value:", round(adf_single_diff$p.value, 4)))

# Plot ACF and PACF Double Differenced, ADF Test
ggtsdisplay(ts_double_diff, main="Double Differenced Data")
adf_double_diff <- adf.test(ts_double_diff)
print(paste("ADF p-value:", round(adf_double_diff$p.value, 4)))

# Individual ACF and PACF Plots - Might be useful for report
ggAcf(ts_double_diff, lag.max=36)+
  theme_minimal(base_size=16)

ggPacf(ts_double_diff, lag.max=36)+
  theme_minimal(base_size=16)

acf(ts_double_diff, lag.max=36)
pacf(ts_double_diff, lag.max=36)

# MODEL FITTING BASELINE =======================================================
# m* = model *
# ARIMA(1,1,1)(0,1,1)[12]
m1 <- Arima(ts_data, order=c(1,1,1), seasonal=c(0,1,1))
print(paste("Model 1 AIC:", round(m1$aic, 2)))
#checkresiduals(m1)

# MODEL FITTING MOVING AVERAGE =================================================
# ARIMA(0,1,1)(0,1,1)[12]
m2 <- Arima(ts_data, order=c(0,1,1), seasonal=c(0,1,1))
print(paste("Model 2 AIC:", round(m2$aic, 2)))
#checkresiduals(m2)

# ARIMA(0,1,2)(0,1,1)[12]
m3 <- Arima(ts_data, order=c(0,1,2), seasonal=c(0,1,1))
print(paste("Model 3 AIC:", round(m3$aic, 2)))
#checkresiduals(m3)

# ARIMA(0,1,3)(0,2,2)[12]
m4 <- Arima(ts_data, order=c(0,1,3), seasonal=c(0,2,2))
print(paste("Model 4 AIC:", round(m4$aic, 2)))
#checkresiduals(m4)

# MODEL FITTING MIXED ==========================================================
# ARIMA(1,1,2)(0,1,1)[12] BEST
m5 <- Arima(ts_data, order=c(1,1,2), seasonal=c(0,1,1))
print(paste("Model 5 AIC:", round(m5$aic, 2)))
#checkresiduals(m5)

# ARIMA(1,1,2)(0,2,2)[12]
m6 <- Arima(ts_data, order=c(1,1,2), seasonal=c(0,2,2))
print(paste("Model 6 AIC:", round(m6$aic, 2)))
#checkresiduals(m6)

# ARIMA(1,1,1)(1,1,0)[12]
m7 <- Arima(ts_data, order=c(1,1,1), seasonal=c(1,1,0))
print(paste("Model 7 AIC:", round(m7$aic, 2)))
#checkresiduals(m7)

# ARIMA(1,1,0)(0,1,1)[12] 
m8 <- Arima(ts_data, order=c(1,1,0), seasonal=c(0,1,1))
print(paste("Model 8 AIC:", round(m8$aic, 2)))
#checkresiduals(m8)

# ARIMA(1,1,3)(0,1,1)[12] 
m9 <- Arima(ts_data, order=c(1,1,3), seasonal=c(0,1,1))
print(paste("Model 9 AIC:", round(m9$aic, 2)))
#checkresiduals(m9)

# MODEL FITTING AUTOREGRESSIVE =================================================
# ARIMA(1,1,0)(1,1,0)[12] 
m10 <- Arima(ts_data, order=c(1,1,0), seasonal=c(1,1,0))
print(paste("Model 10 AIC:", round(m10$aic, 2)))
#checkresiduals(m10)

# SELECTING THE BEST MODEL =====================================================
best_fit <- m5

# LJUNG BOX TEST ===============================================================
# All pass the test
model_list <- list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)
print("Ljung Box Test Results")
for (i in 1:length(model_list)){
  current_m <- model_list[[i]]
  lb_test <- Box.test(residuals(current_m), type = "Ljung-Box")
  print(paste("Model",i, ": p = ", round(lb_test$p.value, 4)))
}

checkresiduals(best_fit)
lb_test_best <- Box.test(residuals(best_fit), type = "Ljung-Box")
print(paste("Best Model - Ljung-Box p-value:", round(lb_test_best$p.value, 4)))

# FORECAST =====================================================================
# Hard coded the model due to issues with using lambda at the same time
final_model <- Arima(ts_data, order=c(1,1,2), seasonal=c(0,1,1), lambda=0)
# Forecast 3 months into 2025
final_forecast <- forecast(final_model, h = 3)
print("FORECAST VALUES Jan-Mar 2025")
print(final_forecast)

# COMPARISON ===================================================================
print("ACTUAL VALUES Jan-Mar 2025")
print(actual_2025)

autoplot(final_forecast, include=6) +
  autolayer(connector_line, color="black") +
  autolayer(actual_2025, series="Actual Data", linewidth=1.5) +
  ggtitle("") +
  ylab("Trade Value (Million of Euro)") +
  labs(caption = "Blue Line: Forecasted Values\nDarker Blue: 80% Confidence Interval\nLighter Blue: 95% Confidence Interval") +
  theme_minimal(base_size=16)

# PERFORMANCE METRICS ==========================================================
performance_metrics <- accuracy(final_forecast, actual_2025)
print("Predictive performance metrics")
print(performance_metrics)
