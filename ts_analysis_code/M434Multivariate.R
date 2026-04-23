library(vars)
library(ggplot2)
library(tseries)
library(forecast)

raw_ts_intex = read.csv("intex_data.csv")
raw_ts_xexp = read.csv("xexp_data.csv")
raw_ts_ximp = read.csv("ximp_data.csv")

# Plotting raw time series =====================================================
raw_cleaned_intex <- subset(raw_ts_intex, geo == "DE")
tv_cleaned_intex <- raw_cleaned_intex$OBS_VALUE
ts_cleaned_intex <- ts(tv_cleaned_intex, start = c(2002, 1), frequency = 12)

autoplot(ts_cleaned_intex) + 
  ggtitle("Intra Export Trade Data 2002 - 2025") + 
  ylab("Trade Value (Million of Euro)") +
  theme_minimal(base_size = 16)

raw_cleaned_xexp <- subset(raw_ts_xexp, geo == "DE")
tv_cleaned_xexp <- raw_cleaned_xexp$OBS_VALUE
ts_cleaned_xexp <- ts(tv_cleaned_xexp, start = c(2002, 1), frequency = 12)

autoplot(ts_cleaned_xexp) + 
  ggtitle("Intra Export Trade Data 2002 - 2025") + 
  ylab("Trade Value (Million of Euro)") +
  theme_minimal(base_size = 16)

raw_cleaned_ximp <- subset(raw_ts_ximp, geo == "DE")
tv_cleaned_ximp <- raw_cleaned_ximp$OBS_VALUE
ts_cleaned_ximp <- ts(tv_cleaned_ximp, start = c(2002, 1), frequency = 12)

autoplot(ts_cleaned_ximp) + 
  ggtitle("Intra Export Trade Data 2002 - 2025") + 
  ylab("Trade Value (Million of Euro)") +
  theme_minimal(base_size = 16)

# Data cleaning ================================================================
clean_intex = raw_ts_intex[raw_ts_intex$geo == "DE", c("TIME_PERIOD", "OBS_VALUE")]
clean_intex$TIME_PERIOD <- as.yearmon(clean_intex$TIME_PERIOD, format = "%Y-%m")
clean_intex = clean_intex[clean_intex$TIME_PERIOD <= as.yearmon("2024-12"), ]
clean_intex <- clean_intex[order(clean_intex$TIME_PERIOD), ]

clean_xexp = raw_ts_xexp[raw_ts_xexp$geo == "DE", c("TIME_PERIOD", "OBS_VALUE")]
clean_xexp$TIME_PERIOD <- as.yearmon(clean_xexp$TIME_PERIOD, format = "%Y-%m")
clean_xexp = clean_xexp[clean_xexp$TIME_PERIOD <= as.yearmon("2024-12"), ]
clean_xexp <- clean_xexp[order(clean_xexp$TIME_PERIOD), ]

clean_ximp = raw_ts_ximp[raw_ts_ximp$geo == "DE", c("TIME_PERIOD", "OBS_VALUE")]
clean_ximp$TIME_PERIOD <- as.yearmon(clean_ximp$TIME_PERIOD, format = "%Y-%m")
clean_ximp = clean_ximp[clean_ximp$TIME_PERIOD <= as.yearmon("2024-12"), ]
clean_ximp <- clean_ximp[order(clean_ximp$TIME_PERIOD), ]

trade_intex<-clean_intex$OBS_VALUE
trade_xexp<-clean_xexp$OBS_VALUE
trade_ximp<-clean_ximp$OBS_VALUE

df_PR<- data.frame(
  intex = trade_intex,
  xexp = trade_xexp,
  ximp = trade_ximp
)
df_PR$month <- rownames(df_PR)

# Actual 2025 Values ===========================================================
actual_2025_intex <- window(ts_cleaned_intex, start = c(2025, 1), end = c(2025, 3))
actual_2025_xexp<- window(ts_cleaned_xexp, start = c(2025, 1), end = c(2025, 3))
actual_2025_ximp <- window(ts_cleaned_ximp, start = c(2025, 1), end = c(2025, 3))

# Initial V1 model
V1 = VAR(cbind(trade_intex, trade_xexp, trade_ximp), type="none")
summary(V1)

V1_pred = predict(V1, n.ahead=3, ci=0.95)
V1_pred

# V(p) model no type
Vp_1 = VAR(cbind(trade_intex, trade_xexp, trade_ximp), type="none", season=12, ic="AIC")
summary(Vp_1)

Vp_1_pred = predict(Vp_1, n.ahead=3, ci=0.95)
Vp_1_pred

# V(p) model with trend
Vp_2 = VAR(cbind(trade_intex, trade_xexp, trade_ximp), type="trend", season=12, ic="AIC")
summary(Vp_2) #eigenvalues less than 1

Vp_2_pred = predict(Vp_2, n.ahead=3, ci=0.95)
Vp_2_pred

df_PR$month <- as.numeric(rownames(df_PR))
n_obs <- nrow(df_PR)

# Plot 1 INTEX =================================================================
fcst_matrix_intex <- Vp_2_pred$fcst$trade_intex 

bridge_intex <- data.frame(
  month = n_obs,
  fcst  = df_PR$intex[n_obs],
  lower = df_PR$intex[n_obs],
  upper = df_PR$intex[n_obs]
)

fcst_rows_intex <- data.frame(
  month = (n_obs + 1):(n_obs + 3),
  fcst  = fcst_matrix_intex[, "fcst"],
  lower = fcst_matrix_intex[, "lower"],
  upper = fcst_matrix_intex[, "upper"]
)

df_intex <- rbind(bridge_intex, fcst_rows_intex)

p1 <- ggplot(df_PR, aes(x=month, y=intex)) +
  geom_line(color="black", linewidth=1) +
  geom_line(data=df_intex, aes(x=month, y=lower), color="steelblue", linetype="dashed", linewidth=0.8) +
  geom_line(data=df_intex, aes(x=month, y=upper), color="steelblue", linetype="dashed", linewidth=0.8) +
  geom_line(data=df_intex, aes(x=month, y=fcst), color="blue", linewidth=1) +
  coord_cartesian(xlim = c(n_obs - 33, n_obs + 3)) + 
  ggtitle("Forecasted Intra Exports") +
  xlab("Month Index") + ylab("Trade Value (Million of Euro)") +
  theme_minimal(base_size = 16)
print(p1)

# Plot 2 XEXP ==================================================================
fcst_matrix_xexp <- PR$fcst$trade_xexp 

bridge_xexp <- data.frame(
  month = n_obs,
  fcst  = df_PR$xexp[n_obs],
  lower = df_PR$xexp[n_obs],
  upper = df_PR$xexp[n_obs]
)

fcst_rows_xexp <- data.frame(
  month = (n_obs + 1):(n_obs + 3),
  fcst  = fcst_matrix_xexp[, "fcst"],
  lower = fcst_matrix_xexp[, "lower"],
  upper = fcst_matrix_xexp[, "upper"]
)

df_xexp <- rbind(bridge_xexp, fcst_rows_xexp)

p2 <- ggplot(df_PR, aes(x=month, y=xexp)) +
  geom_line(color="black", linewidth=1) +
  geom_line(data=df_xexp, aes(x=month, y=lower), color="steelblue", linetype="dashed", linewidth=0.8) +
  geom_line(data=df_xexp, aes(x=month, y=upper), color="steelblue", linetype="dashed", linewidth=0.8) +
  geom_line(data=df_xexp, aes(x=month, y=fcst), color="blue", linewidth=1) +
  coord_cartesian(xlim = c(n_obs - 33, n_obs + 3)) + 
  ggtitle("Forecasted Extra Exports") +
  xlab("Month Index") + ylab("Trade Value (Million of Euro)") +
  theme_minimal(base_size = 16)
print(p2)

# Plot 3 XIMP ==================================================================
fcst_matrix_ximp <- PR$fcst$trade_ximp 

bridge_ximp <- data.frame(
  month = n_obs,
  fcst  = df_PR$ximp[n_obs],
  lower = df_PR$ximp[n_obs],
  upper = df_PR$ximp[n_obs]
)

fcst_rows_ximp <- data.frame(
  month = (n_obs + 1):(n_obs + 3),
  fcst  = fcst_matrix_ximp[, "fcst"],
  lower = fcst_matrix_ximp[, "lower"],
  upper = fcst_matrix_ximp[, "upper"]
)

df_ximp <- rbind(bridge_ximp, fcst_rows_ximp)

p3 <- ggplot(df_PR, aes(x=month, y=ximp)) +
  geom_line(color="black", linewidth=1) +
  geom_line(data=df_ximp, aes(x=month, y=lower), color="steelblue", linetype="dashed", linewidth=0.8) +
  geom_line(data=df_ximp, aes(x=month, y=upper), color="steelblue", linetype="dashed", linewidth=0.8) +
  geom_line(data=df_ximp, aes(x=month, y=fcst), color="blue", linewidth=1) +
  coord_cartesian(xlim = c(n_obs - 33, n_obs + 3)) + 
  ggtitle("Forecasted Extra Imports") +
  xlab("Month Index") + ylab("Trade Value (Million of Euro)") +
  theme_minimal(base_size = 16)
print(p3)

# PERFORMANCE METRICS ==========================================================
print(fcst_intex)
print(actual_2025_intex)
print(fcst_xexp)
print(actual_2025_xexp)
print(fcst_ximp)
print(actual_2025_ximp)

fcst_intex <- fcst_matrix_intex[, "fcst"]
performance_metrics_intex <- accuracy(fcst_intex, actual_2025_intex)
print("Predictive performance metrics intra exports")
print(performance_metrics_intex)

fcst_xexp <- fcst_matrix_xexp[, "fcst"]
performance_metrics_xexp <- accuracy(fcst_xexp, actual_2025_xexp)
print("Predictive performance metrics extra exports")
print(performance_metrics_xexp)

fcst_ximp <- fcst_matrix_ximp[, "fcst"]
performance_metrics_ximp <- accuracy(fcst_ximp, actual_2025_ximp)
print("Predictive performance metrics extra imports")
print(performance_metrics_ximp)


