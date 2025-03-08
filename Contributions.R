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
df <- read.csv("data/contributions-global-temp-change.csv")
View(df)
attach(df)
head(df)
GW = Share.of.contribution.to.global.warming
USA = GW[38926:39098]
China = GW[7613:7785]
European_Union = GW[12976:13148]
Russia = GW[30622:30794]
Brazil = GW[5364:5536]
India = GW[17128:17300]
United_kingdom = GW[38753:38925]
South_America = GW[34428:34600]
Africa = GW[174:346]
Date = Year


USA_ts <- ts(USA, frequency = 1, start = c(1850, 1))
China_ts <- ts(China, frequency = 1, start = c(1850, 1))
European_Union_ts <- ts(European_Union, frequency = 1, start = c(1850, 1))
Russia_ts <- ts(Russia, frequency = 1, start = c(1850, 1))
Brazil_ts <- ts(Brazil, frequency = 1, start = c(1850, 1))
India_ts <- ts(India, frequency = 1, start = c(1850, 1))
United_kingdom_ts <- ts(United_kingdom, frequency = 1, start = c(1850, 1))
South_America_ts <- ts(South_America, frequency = 1, start = c(1850, 1))
Africa_ts <- ts(Africa, frequency = 1, start = c(1850, 1))


countries_data <- data.frame(
  Year = c(time(USA_ts), time(China_ts), time(European_Union_ts), time(Russia_ts),
           time(Brazil_ts), time(India_ts), time(United_kingdom_ts), time(South_America_ts), time(Africa_ts)),
  Contribution = c(as.numeric(USA_ts), as.numeric(China_ts), as.numeric(European_Union_ts),
                   as.numeric(Russia_ts), as.numeric(Brazil_ts), as.numeric(India_ts),
                   as.numeric(United_kingdom_ts), as.numeric(South_America_ts), as.numeric(Africa_ts)),
  Country = factor(rep(c("USA", "China", "European Union", "Russia", "Brazil",
                         "India", "United Kingdom", "South America", "Africa"),
                       times = c(length(USA_ts), length(China_ts), length(European_Union_ts),
                                 length(Russia_ts), length(Brazil_ts), length(India_ts),
                                 length(United_kingdom_ts), length(South_America_ts), length(Africa_ts))))
)

ggplot(countries_data, aes(x = Year, y = Contribution, color = Country)) +
  geom_line(size = 1) +
  ggtitle("Yearly Contribution to Global Warming by Country") +
  xlab("Year") +
  ylab("Yearly Share of Contribution to Global Warming") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("USA" = "blue", "China" = "red", "European Union" = "green",
                                "Russia" = "purple", "Brazil" = "orange", "India" = "brown",
                                "United Kingdom" = "pink", "South America" = "cyan", "Africa" = "darkgreen"))

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
European_Union_result <- perform_regression(European_Union, "European Union")
Russia_result <- perform_regression(Russia, "Russia")
Brazil_result <- perform_regression(Brazil, "Brazil")
India_result <- perform_regression(India, "India")
United_kingdom_result <- perform_regression(United_kingdom, "United Kingdom")
South_America_result <- perform_regression(South_America, "South America")
Africa_result <- perform_regression(Africa, "Africa")
####################################################################################################################
## TSLM model
####################################################################################################################
perform_tslm <- function(data, country_name) {
  ts_data <- ts(data, frequency = 1, start = c(1850, 1))

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

# Apply the TSLM function to each country
USA_result <- perform_tslm(USA, "USA")
China_result <- perform_tslm(China, "China")
European_Union_result <- perform_tslm(European_Union, "European Union")
Russia_result <- perform_tslm(Russia, "Russia")
Brazil_result <- perform_tslm(Brazil, "Brazil")
India_result <- perform_tslm(India, "India")
United_kingdom_result <- perform_tslm(United_kingdom, "United Kingdom")
South_America_result <- perform_tslm(South_America, "South America")
Africa_result <- perform_tslm(Africa, "Africa")
######################################################################################################################
## Exponential Smoothing Methods
######################################################################################################################
perform_ses <- function(data, country_name, alpha_values = c(0.2, 0.6), forecast_horizon = 5) {
  ts_data <- ts(data, frequency = 1, start = c(1850, 1)) # Adjust frequency and start year as needed

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
European_Union_ses <- perform_ses(European_Union, "European Union")
Russia_ses <- perform_ses(Russia, "Russia")
Brazil_ses <- perform_ses(Brazil, "Brazil")
India_ses <- perform_ses(India, "India")
United_kingdom_ses <- perform_ses(United_kingdom, "United Kingdom")
South_America_ses <- perform_ses(South_America, "South America")
Africa_ses <- perform_ses(Africa, "Africa")
##########################################################################################################################
## Trend methods (Holt method)
##########################################################################################################################
perform_holt <- function(data, country_name, forecast_horizon = 5) {
  ts_data <- ts(data, frequency = 1, start = c(1850, 1))
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
European_Union_holt <- perform_holt(European_Union, "European Union")
Russia_holt <- perform_holt(Russia, "Russia")
Brazil_holt <- perform_holt(Brazil, "Brazil")
India_holt <- perform_holt(India, "India")
United_kingdom_holt <- perform_holt(United_kingdom, "United Kingdom")
South_America_holt <- perform_holt(South_America, "South America")
Africa_holt <- perform_holt(Africa, "Africa")

##########################################################################################################################
## Arima Models
##########################################################################################################################
perform_arima <- function(data, country_name, forecast_horizon = 5, lag=3) {
  ts_data <- ts(data, frequency = 1, start = c(1850, 1))
  arima_model <- auto.arima(ts_data)
  fc_arima <- forecast(arima_model, h = forecast_horizon)
  cat("5 Years Forecast for  ", country_name, ":\n")
  print(fc_arima)
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
  box_pierce <- Box.test(residuals_arima, type = "Box-Pierce", lag = lag)
  return(list(arima_model = arima_model, fc_arima = fc_arima, accuracy_arima = accuracy_arima, box_pierce_result = box_pierce))
}

USA_arima <- perform_arima(USA, "USA")
China_arima <- perform_arima(China, "China")
European_Union_arima <- perform_arima(European_Union, "European Union")
Russia_arima <- perform_arima(Russia, "Russia", lag = 1)
Brazil_arima <- perform_arima(Brazil, "Brazil", lag = 2)
India_arima <- perform_arima(India, "India")
United_kingdom_arima <- perform_arima(United_kingdom, "United Kingdom")
South_America_arima <- perform_arima(South_America, "South America", lag = 1)
Africa_arima <- perform_arima(Africa, "Africa")

##############################################################################################################################
## Arima model plot
##############################################################################################################################
combine_forecast_data <- function(actual_data, forecasts, country_names, start_year_actual = 1850) {
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
  actual_data = list(USA, China, European_Union, Russia, Brazil, India, United_kingdom, South_America, Africa),
  forecasts = list(USA_arima$fc_arima, China_arima$fc_arima, European_Union_arima$fc_arima,
                   Russia_arima$fc_arima, Brazil_arima$fc_arima, India_arima$fc_arima,
                   United_kingdom_arima$fc_arima, South_America_arima$fc_arima, Africa_arima$fc_arima),
  country_names = c("USA", "China", "European Union", "Russia", "Brazil", "India",
                    "United Kingdom", "South America", "Africa")
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
