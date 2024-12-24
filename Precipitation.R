library(readxl)
library(lmtest)
library(forecast)
library(DIMORA)
library(fpp2)
library(dplyr)
library(car)
#############################################################################################################################
# Data presentation
#############################################################################################################################
df <- read.csv("data/average-precipitation-per-year.csv")
View(df)
attach(df)
head(df)
AP = Annual.precipitation
USA = AP[15457:15540]
China = AP[3109:3192]
Russia = AP[12265:12348]
Brazil = AP[2017:2100]
India = AP[7057:7140]
United_kingdom = AP[15373:15456]
Date = Year


USA_ts <- ts(USA, frequency = 1, start = c(1940, 1))
China_ts <- ts(China, frequency = 1, start = c(1940, 1))
Russia_ts <- ts(Russia, frequency = 1, start = c(1940, 1))
Brazil_ts <- ts(Brazil, frequency = 1, start = c(1940, 1))
India_ts <- ts(India, frequency = 1, start = c(1940, 1))
United_kingdom_ts <- ts(United_kingdom, frequency = 1, start = c(1940, 1))


countries_data <- data.frame(
  Year = c(time(USA_ts), time(China_ts), time(Russia_ts),
           time(Brazil_ts), time(India_ts), time(United_kingdom_ts)),
  Contribution = c(as.numeric(USA_ts), as.numeric(China_ts),
                   as.numeric(Russia_ts), as.numeric(Brazil_ts), as.numeric(India_ts),
                   as.numeric(United_kingdom_ts)),
  Country = factor(rep(c("USA", "China", "Russia", "Brazil",
                         "India", "United Kingdom"),
                       times = c(length(USA_ts), length(China_ts),
                                 length(Russia_ts), length(Brazil_ts), length(India_ts),
                                 length(United_kingdom_ts))))
)

ggplot(countries_data, aes(x = Year, y = Contribution, color = Country)) +
  geom_line(size = 1) +
  ggtitle("Annual Precipitation for Top CO2 and Methane Producing Countries") +
  xlab("Year") +
  ylab("Annual Precipitation") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c(
    "USA" = "blue",
    "China" = "red",
    "Russia" = "purple",
    "Brazil" = "orange",
    "India" = "brown",
    "United Kingdom" = "pink"
  ))
###################################################################################################################
## Linear Regression
###################################################################################################################

perform_regression <- function(data, country_name) {
  tt <- 1:NROW(data)

  model <- lm(data ~ tt)

  cat("Linear Regression Summary for", country_name, "\n")
  print(summary(model))
  print(anova(model))

  fitted_values <- fitted(model)
  residuals_values <- residuals(model)
  mae <- mean(abs(residuals_values))
  rmse <- sqrt(mean(residuals_values^2))
  cat("Accuracy Metrics for", country_name, ":\n")
  cat("MAE:", mae, "\nRMSE:", rmse, "\n\n")

  dw_result <- dwtest(model)
  cat("Durbin-Watson Test for", country_name, ":\n")
  print(dw_result)

  par(mfrow = c(2, 2))
  plot(tt, data, xlab = "Year", ylab = "Contribution to Global Warming",
       main = paste("Linear Regression for", country_name))
  abline(model, col = "blue", lwd = 2)

  plot(residuals_values, xlab = "Year", ylab = "Residuals",
       main = paste("Residuals for", country_name))

  Acf(residuals_values, main = paste("ACF of Residuals for", country_name), lag.max = 180)
  pacf(residuals_values, main = paste("PACF of Residuals for", country_name), lag.max = 180)
  par(mfrow = c(1, 1))

  return(list(model = model, residuals = residuals_values))
}


USA_result <- perform_regression(USA, "USA")
China_result <- perform_regression(China, "China")
Russia_result <- perform_regression(Russia, "Russia")
Brazil_result <- perform_regression(Brazil, "Brazil")
India_result <- perform_regression(India, "India")
United_kingdom_result <- perform_regression(United_kingdom, "United Kingdom")
####################################################################################################################
## TSLM model
####################################################################################################################
perform_tslm <- function(data, country_name) {
  ts_data <- ts(data, frequency = 1, start = c(1940, 1))

  TSLM_model <- tslm(ts_data ~ trend)

  cat("TSLM Summary for", country_name, "\n")
  print(summary(TSLM_model))
  tslm_accuracy <- accuracy(TSLM_model)
  cat("Accuracy Metrics for", country_name, ":\n")
  print(tslm_accuracy)

  par(mfrow = c(3, 2))
  plot(ts_data, main = paste("TSLM Model for", country_name),
       xlab = "Year", ylab = "Contribution to Global Warming")
  lines(fitted(TSLM_model), col = "red", lwd = 2)
  legend("topleft", legend = c("Original Data", "Fitted Model"), col = c("black", "red"), lty = 1, lwd = 2)

  TSLM_res <- residuals(TSLM_model)
  plot(TSLM_res, main = paste("Residuals for", country_name), xlab = "Year", ylab = "Residuals")
  Acf(TSLM_res, main = paste("ACF of Residuals for", country_name), lag.max = 180)
  pacf(TSLM_res, main = paste("PACF of Residuals for", country_name), lag.max = 180)

  dw_result <- dwtest(TSLM_model)
  cat("Durbin-Watson Test for", country_name, ":\n")
  print(dw_result)

  TSLM_fore <- forecast(TSLM_model, h = 5)
  plot(TSLM_fore, main = paste("Forecast for", country_name), xlab = "Year", ylab = "Contribution to Global Warming")
  par(mfrow = c(1, 1))

  return(list(model = TSLM_model, accuracy = tslm_accuracy, residuals = TSLM_res))
}

USA_result <- perform_tslm(USA, "USA")
China_result <- perform_tslm(China, "China")
Russia_result <- perform_tslm(Russia, "Russia")
Brazil_result <- perform_tslm(Brazil, "Brazil")
India_result <- perform_tslm(India, "India")
United_kingdom_result <- perform_tslm(United_kingdom, "United Kingdom")
######################################################################################################################
## Exponential Smoothing Methods
######################################################################################################################
perform_ses <- function(data, country_name, alpha_values = c(0.2, 0.6), forecast_horizon = 5) {
  ts_data <- ts(data, frequency = 1, start = c(1940, 1)) # Adjust frequency and start year as needed

  ses_fit1 <- ses(ts_data, alpha = alpha_values[1], initial = "simple", h = forecast_horizon)
  ses_fit2 <- ses(ts_data, alpha = alpha_values[2], initial = "simple", h = forecast_horizon)
  ses_fit3 <- ses(ts_data, h = forecast_horizon)

  par(mfrow = c(3, 2))
  plot(ts_data, main = paste("SES Model 1 for", country_name),
       ylab = "Contribution to Global Warming", xlab = "Time (Month of Year)")
  lines(fitted(ses_fit1), col = "blue", lwd = 2)
  legend("topleft", legend = c("Original Data", "Fitted Model 1"), col = c("black", "blue"), lty = 1, lwd = 2)

  plot(ts_data, main = paste("SES Model 2 for", country_name),
       ylab = "Contribution to Global Warming", xlab = "Year")
  lines(fitted(ses_fit2), col = "red", lwd = 2)
  legend("topleft", legend = c("Original Data", "Fitted Model 2"), col = c("black", "red"), lty = 1, lwd = 2)

  plot(ts_data, main = paste("SES Model 3 for", country_name, "(Automatic Alpha)"),
       ylab = "Contribution to Global Warming", xlab = "Year")
  lines(fitted(ses_fit3), col = "green", lwd = 2)
  legend("topleft", legend = c("Original Data", "Fitted Model 3"), col = c("black", "green"), lty = 1, lwd = 2)


  ses_forecast <- forecast(ses_fit3, h = forecast_horizon)

  plot(ses_forecast, main = paste("Forecast for", country_name),
       xlab = "Year", ylab = "Contribution to Global Warming",
       col = "blue", lwd = 2)
  lines(fitted(ses_fit3), col = "red", lwd = 2)
  legend("topleft", legend = c("Forecast", "Fitted"), col = c("blue", "red"), lty = 1, lwd = 2)

  accuracy_fit3 <- accuracy(ses_fit3)
  cat("Accuracy Metrics for SES (Fit 3) -", country_name, ":\n")
  print(accuracy_fit3)

  residuals_fit3 <- residuals(ses_fit3)
  cat("Residual Summary for SES (Fit 3) -", country_name, ":\n")
  print(summary(residuals_fit3))

  Acf(residuals_fit3, main = paste("ACF of Residuals for SES (Fit 3) -", country_name), lag.max = 180)
  pacf(residuals_fit3, main = paste("PACF of Residuals for SES (Fit 3) -", country_name), lag.max = 180)
  par(mfrow = c(1, 1))

  return(list(fit1 = ses_fit1, fit2 = ses_fit2, fit3 = ses_fit3, forecast = ses_forecast, accuracy = accuracy_fit3))
}

USA_ses <- perform_ses(USA, "USA")
China_ses <- perform_ses(China, "China")
Russia_ses <- perform_ses(Russia, "Russia")
Brazil_ses <- perform_ses(Brazil, "Brazil")
India_ses <- perform_ses(India, "India")
United_kingdom_ses <- perform_ses(United_kingdom, "United Kingdom")
##########################################################################################################################
## Trend methods (Holt method)
##########################################################################################################################
perform_holt <- function(data, country_name, forecast_horizon = 5) {
  ts_data <- ts(data, frequency = 1, start = c(1940, 1))
  fc_Holt <- holt(ts_data, h = forecast_horizon)
  cat("5 Years Forecast for  ", country_name, ":\n")
  print(fc_Holt)
  par(mfrow = c(3, 1))
  plot(fc_Holt, main = paste("Holt's Method for", country_name),
       xlab = "Year", ylab = "Contribution to Global Warming", col = "blue", lwd = 2)
  lines(fitted(fc_Holt), col = "red", lwd = 2)
  legend("topleft", legend = c("Forecast", "Fitted"), col = c("blue", "red"), lty = 1, lwd = 2)
  accuracy_fc_Holt <- accuracy(fc_Holt)
  cat("Accuracy Metrics for Holt's Method -", country_name, ":\n")
  print(accuracy_fc_Holt)
  residuals_Holt <- residuals(fc_Holt)
  cat("Residual Summary for Holt's Method -", country_name, ":\n")
  print(summary(residuals_Holt))
  Acf(residuals_Holt, main = paste("ACF of Residuals for Holt's Method -", country_name), lag.max = 180)
  pacf(residuals_Holt, main = paste("PACF of Residuals for Holt's Method -", country_name), lag.max = 180)
  par(mfrow = c(1, 1))
  return(list(fc_Holt = fc_Holt, accuracy_fc_Holt = accuracy_fc_Holt))
}

USA_holt <- perform_holt(USA, "USA")
China_holt <- perform_holt(China, "China")
Russia_holt <- perform_holt(Russia, "Russia")
Brazil_holt <- perform_holt(Brazil, "Brazil")
India_holt <- perform_holt(India, "India")
United_kingdom_holt <- perform_holt(United_kingdom, "United Kingdom")

##########################################################################################################################
## Arima Models
##########################################################################################################################
# perform_arima <- function(data, country_name, forecast_horizon = 5) {
#   ts_data <- ts(data, frequency = 1, start = c(1940, 1))
#   ts_data <- diff(ts_data)
#   ts_data <- diff(ts_data, lag = 1)
#   arima_model <- auto.arima(ts_data)
#   fc_arima <- forecast(arima_model, h = forecast_horizon)
#   cat("5 Years Forecast for  ", country_name, ":\n")
#   print(fc_arima)
#   par(mfrow = c(3,1))
#   plot(fc_arima, main = paste("ARIMA Model for", country_name),
#        xlab = "Year", ylab = "Contribution to Global Warming", col = "blue", lwd = 2)
#   lines(fitted(arima_model), col = "red", lwd = 2)
#   legend("topleft", legend = c("Forecast", "Fitted"), col = c("blue", "red"), lty = 1, lwd = 2)
#   accuracy_arima <- accuracy(fc_arima)
#   cat("Accuracy Metrics for ARIMA Model -", country_name, ":\n")
#   print(accuracy_arima)
#   residuals_arima <- residuals(arima_model)
#   cat("Residual Summary for ARIMA Model -", country_name, ":\n")
#   print(summary(residuals_arima))
#   Acf(residuals_arima, main = paste("ACF of Residuals for ARIMA Model -", country_name),lag.max = 180)
#   pacf(residuals_arima, main = paste("PACF of Residuals for ARIMA Model -", country_name),lag.max = 180)
#   par(mfrow = c(1, 1))
#   return(list(arima_model = arima_model, fc_arima = fc_arima, accuracy_arima = accuracy_arima))
# }

perform_arima <- function(data, country_name, forecast_horizon = 5) {
  ts_data <- ts(data, frequency = 1, start = c(1940, 1))
  # Differencing to make data stationary
  ts_data_diff <- diff(ts_data)
  ts_data_diff <- diff(ts_data_diff, lag = 1)
  # Fit ARIMA on differenced data
  arima_model <- auto.arima(ts_data_diff)
  fc_arima <- forecast(arima_model, h = forecast_horizon)
  cat("5 Years Forecast for  ", country_name, ":\n")
  print(fc_arima)
  # Inverse differencing to return forecasts to original scale
  fc_arima$mean <- cumsum(c(ts_data[length(ts_data)], fc_arima$mean)) # Adding last actual value for correct forecast alignment

  par(mfrow = c(3, 1))
  plot(fc_arima, main = paste("ARIMA Model for", country_name),
       xlab = "Year", ylab = "Contribution to Global Warming", col = "blue", lwd = 2)
  lines(fitted(arima_model), col = "red", lwd = 2)
  legend("topleft", legend = c("Forecast", "Fitted"), col = c("blue", "red"), lty = 1, lwd = 2)

  accuracy_arima <- accuracy(fc_arima)
  cat("Accuracy Metrics for ARIMA Model -", country_name, ":\n")
  print(accuracy_arima)

  residuals_arima <- residuals(arima_model)
  cat("Residual Summary for ARIMA Model -", country_name, ":\n")
  print(summary(residuals_arima))
  Acf(residuals_arima, main = paste("ACF of Residuals for ARIMA Model -", country_name), lag.max = 180)
  pacf(residuals_arima, main = paste("PACF of Residuals for ARIMA Model -", country_name), lag.max = 180)

  par(mfrow = c(1, 1))
  return(list(arima_model = arima_model, fc_arima = fc_arima, accuracy_arima = accuracy_arima))
}

USA_arima <- perform_arima(USA_ts, "USA")
China_arima <- perform_arima(China_ts, "China")
Russia_arima <- perform_arima(Russia_ts, "Russia")
Brazil_arima <- perform_arima(Brazil_ts, "Brazil")
India_arima <- perform_arima(India_ts, "India")
United_kingdom_arima <- perform_arima(United_kingdom_ts, "United Kingdom")
##############################################################################################################################
## Arima model plot
##############################################################################################################################
combine_forecast_data <- function(actual_data, forecasts, country_names, start_year_actual = 1940) {
  combined_data <- data.frame()

  for (i in seq_along(forecasts)) {
    actual_years <- seq(start_year_actual, by = 1, length.out = length(actual_data[[i]]))
    actual_df <- data.frame(
      Year = actual_years,
      Contribution = actual_data[[i]],
      Country = country_names[i],
      Type = "Actual"
    )

    # Extract the point forecast (e.g., mean or specific value)
    point_forecast <- if (is.list(forecasts[[i]])) forecasts[[i]]$mean else forecasts[[i]]

    forecast_years <- seq(max(actual_years) + 1, by = 1, length.out = length(point_forecast))
    forecast_df <- data.frame(
      Year = forecast_years,
      Contribution = as.numeric(point_forecast),
      Country = country_names[i],
      Type = "Forecast"
    )

    combined_data <- rbind(combined_data, actual_df, forecast_df)
  }

  return(combined_data)
}

# Data for ARIMA forecasts
combined_arima_data <- combine_forecast_data(
  actual_data = list(USA, China, Russia, Brazil, India, United_kingdom),
  forecasts = list(USA_arima$fc_arima, China_arima$fc_arima,
                   Russia_arima$fc_arima, Brazil_arima$fc_arima, India_arima$fc_arima,
                   United_kingdom_arima$fc_arima),
  country_names = c("USA", "China", "Russia", "Brazil", "India",
                    "United Kingdom")
)

# Plot for ARIMA forecasts
ggplot(combined_arima_data, aes(x = Year, y = Contribution, color = Country, linetype = Type)) +
  geom_line(size = 1) +
  ggtitle("Yearly Contribution to Global Warming with ARIMA Point Forecasts") +
  xlab("Year") +
  ylab("Yearly Share of Contribution to Global Warming") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  scale_linetype_manual(values = c("Actual" = "solid", "Forecast" = "dashed")) +
  theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10))

#############################################################################################################################
## Holt Model plot
#############################################################################################################################
combined_holt_data <- combine_forecast_data(
  actual_data = list(USA, China, European_Union, Russia, Brazil, India, United_kingdom, South_America, Africa),
  forecasts = list(USA_holt$fc_Holt, China_holt$fc_Holt, European_Union_holt$fc_Holt,
                   Russia_holt$fc_Holt, Brazil_holt$fc_Holt, India_holt$fc_Holt,
                   United_kingdom_holt$fc_Holt, South_America_holt$fc_Holt, Africa_holt$fc_Holt),
  country_names = c("USA", "China", "European Union", "Russia", "Brazil", "India",
                    "United Kingdom", "South America", "Africa")
)

# Plot for Holt forecasts
ggplot(combined_holt_data, aes(x = Year, y = Contribution, color = Country, linetype = Type)) +
  geom_line(size = 1) +
  ggtitle("Yearly Contribution to Global Warming with Holt's Point Forecasts") +
  xlab("Year") +
  ylab("Yearly Share of Contribution to Global Warming") +
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +
  scale_linetype_manual(values = c("Actual" = "solid", "Forecast" = "dashed")) +
  theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10))