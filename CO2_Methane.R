library(readxl)
library(lmtest) 
library(forecast)
library(DIMORA)
library(fpp2)
library(dplyr)
library(car)
library(gridExtra) 
############################################################################################################
# Data presentation
############################################################################################################
df <- read.csv("data/global-co2-concentration.csv")
View(df)
attach(df)
head(df)
CO2 =  Monthly.concentration.of.atmospheric.carbon.dioxide
Date = Day
tsdisplay(CO2, lag.max = 600)
CO2_ts <- ts(CO2, frequency = 12, start = c(1979, 1))
seasonplot(CO2_ts, ylab = "Monthly concentration of Atmospheric Carbon Dioxide", xlab = "Month of Year",
           main = "Seasonal Plot: concentration of Atmospheric Carbon Dioxide", 
           year.labels = TRUE, year.labels.left = TRUE, 
           col = 1:40, pch = 19)


decomposition <- stl(CO2_ts, s.window = 'periodic')


components <- data.frame(
  Date = as.Date(time(CO2_ts)), # Convert time series to dates
  Observed = as.numeric(CO2_ts), # Original data
  Seasonal = decomposition$time.series[, "seasonal"],
  Trend = decomposition$time.series[, "trend"],
  Remainder = decomposition$time.series[, "remainder"]
)

p1 <- ggplot(components, aes(x = Date, y = Observed)) +
  geom_line(color = "blue") +
  labs(title = "Observed CO2", x = "", y = "CO2") +
  theme_minimal()

p2 <- ggplot(components, aes(x = Date, y = Seasonal)) +
  geom_line(color = "orange") +
  labs(title = "Seasonal Component", x = "", y = "Seasonal") +
  theme_minimal()

p3 <- ggplot(components, aes(x = Date, y = Trend)) +
  geom_line(color = "green") +
  labs(title = "Trend Component", x = "", y = "Trend") +
  theme_minimal()

p4 <- ggplot(components, aes(x = Date, y = Remainder)) +
  geom_line(color = "red") +
  labs(title = "Remainder", x = "Date", y = "Remainder") +
  theme_minimal()

grid.arrange(p1, p2, p3, p4, ncol = 1)




################################################################################################################
# Linear Regression
################################################################################################################
tt<- 1:NROW(df)
model1 <- lm(CO2~ tt)
summary(model1)
anova(model1)
accuracy(model1)

##plot of the model
plot(tt, CO2, xlab="Time(Month of Year)", ylab="Monthly concentration of Atmospheric Carbon Dioxide")
abline(model1, col=6, lwd = 2)

dwtest(model1)

resmodel1<- residuals(model1)
plot(resmodel1,xlab="Time(Month of Year)", ylab="residuals" )
Acf(resmodel1, lag.max = 22*12)
pacf(resmodel1, lag.max = 22*12)
###################################################################################################################
# TSLM model  
###################################################################################################################
ts.plot(CO2_ts, type="o")

## we fit a linear model with the tslm function
TSLM_model<- tslm(CO2_ts~ trend + season, data = df)

###obviously it gives the same results of the first model
summary(TSLM_model)
accuracy(TSLM_model)

par(mfrow = c(1, 2)) 
par(cex.lab = 0.7, cex.axis = 0.7, cex.main = 0.6) 

plot(CO2_ts, 
     xlab = "Time (Month of Year)", 
     ylab = "Monthly concentration of Atmospheric Carbon Dioxide",
     main = "Time Series with Fitted Values")
lines(fitted(TSLM_model), col = 2)


TSLM_fore <- forecast(TSLM_model, h = 12)
plot(TSLM_fore, 
     main = "Forecast for the Next 12 Months for Atmospheric Carbon Dioxide", 
     xlab = "Time (Month of Year)", 
     ylab = "Monthly concentration of Atmospheric Carbon Dioxide")

par(mfrow = c(1, 1))

TSLM_res<- residuals(TSLM_model)
par(mfrow = c(1,3))
par(cex.lab = 0.9, cex.axis = 0.9, cex.main = 0.9) 
plot(TSLM_res,xlab="Time(Month of Year)", ylab="residuals", main = "TSLM Model Residuals Plot")
Acf(TSLM_res, lag.max = 22*12,xlab="Time(Month of Year)", ylab="ACF for Atmospheric Carbon Dioxide" , main = "ACF for TSLM Model")
pacf(TSLM_res, lag.max = 22*12,xlab="Time(Month of Year)", ylab="PACF for Atmospheric Carbon Dioxide",main = "PACF for TSLM Model")
par(mfrow = c(1,1))

dwtest(TSLM_model)


###################################################################################################################
##Exponential Smoothing methods
###################################################################################################################
##1.Simple exponential smoothing
autoplot(CO2_ts)+ylab("Monthly concentration of Atmospheric Carbon Dioxide")+xlab("Time(Month of Year)")

fit1<- ses(CO2_ts, alpha=0.2, initial="simple", h=5)
fit2<- ses(CO2_ts, alpha=0.6, initial="simple", h=5)
fit3<- ses(CO2_ts, h=5)


par(mfrow = c(3,1))

plot(CO2_ts, ylab = "Monthly concentration of Atmospheric Carbon Dioxide", 
     xlab = "Time (Month of Year)", main = "Original Time Series with Fit 1")
lines(fitted(fit1), col = "blue", type = "o", lwd = 2)

plot(CO2_ts, ylab = "Monthly concentration of Atmospheric Carbon Dioxide", 
     xlab = "Time (Month of Year)", main = "Original Time Series with Fit 2")
lines(fitted(fit2), col = "red", type = "o", lwd = 2)


plot(CO2_ts, ylab = "Monthly concentration of Atmospheric Carbon Dioxide", 
     xlab = "Time (Month of Year)", main = "Original Time Series with Fit 3")
lines(fitted(fit3), col = "green", type = "o", lwd = 2)
par(mfrow = c(1,1))


fc<- ses(CO2_ts, h=12)
accuracy(fc)

summary(fc)

autoplot(fc)+
  autolayer(fitted(fc), series="Fitted")+ylab("Monthly concentration of Atmospheric Carbon Dioxide")+xlab("Time (Month of Year)")



##2.Trend methods (Holt method)
fc_Holt<- holt(CO2_ts, h=12)
fc2_Holt<- holt(CO2_ts, damped=T, phi=0.9, h=12)

autoplot(CO2_ts)+
  autolayer(fc_Holt, series="Holt's method", PI=F)+
  autolayer(fc2_Holt, series="Damped Holt's method", PI=F)

###3.Trend and seasonality methods (Holt-Winters method)
autoplot(CO2_ts)

fit1<- hw(CO2_ts, seasonal="additive")
fit2<- hw(CO2_ts, seasonal="multiplicative")

autoplot(CO2_ts, series = "Original Time Series") +
  autolayer(fit1, series = "Additive Holt-Winters", PI = FALSE) +
  autolayer(fit2, series = "Multiplicative Holt-Winters", PI = FALSE) +
  ggtitle("Holt-Winters Forecasts with Different Seasonal Models") +
  xlab("Time (Month of Year)") +
  ylab("Monthly Concentration of Atmospheric Carbon Dioxide") +
  theme_minimal() +
  scale_color_manual(
    values = c("Original Time Series" = "black",
               "Additive Holt-Winters" = "blue", 
               "Multiplicative Holt-Winters" = "red")
  ) +
  theme(
    plot.title = element_text(size = 8),        
    axis.title.x = element_text(size = 8),      
    axis.title.y = element_text(size = 8),      
    axis.text = element_text(size = 8),
    legend.text = element_text(size = 8),      
    legend.title = element_text(size = 8) 
  )

  
# Evaluate performance of fit1 and fit2
accuracy_fit1 <- accuracy(fit1)
accuracy_fit2 <- accuracy(fit2)

par(mfrow = c(2,2))
par(cex.lab = 0.7, cex.axis = 0.7, cex.main = 0.7) 
res_fit1 <- residuals(fit1)
Acf(res_fit1, main = "ACF of Residuals (Additive)", lag.max = 22*12)
pacf(res_fit1, main = "PACF of Residuals (Additive)",lag.max = 22*12)

res_fit2 <- residuals(fit2)
Acf(res_fit2, main = "ACF of Residuals (Multiplicative)",lag.max = 22*12)
pacf(res_fit2, main = "PACF of Residuals (Multiplicative)",lag.max = 22*12)
par(mfrow = c(1,1))

########################################################################################################
# SARIMA Models
########################################################################################################
diff1 <- diff(CO2_ts)

# Seasonal differencing to remove the seasonal pattern
diff_seasonal <- diff(diff1, lag = 12)

# Display ACF and PACF after differencing
tsdisplay(diff_seasonal, lag.max = 22*12)


# first arima model
sarima_model <- Arima(CO2_ts, order = c(1, 1, 1), seasonal = list(order = c(0, 1, 0), period = 12))

fit1<- fitted(sarima_model)

par(mfrow = c(1, 2))
par(cex.lab = 0.7, cex.axis = 0.7, cex.main = 0.7) 
plot(CO2_ts, 
     xlab = "Time (Month of Year)", 
     ylab = "Monthly concentration of Atmospheric Carbon Dioxide",
     main = "Time Series with Fitted Values")
lines(fitted(sarima_model), col = 2)
legend("topleft", 
       legend = c("Observed", "Fitted (SARIMA)"), 
       col = c("black", "red"), 
       lty = c(1, 1), 
       bty = "n", 
       cex = 0.8)


sarima_fore <- forecast(sarima_model, h = 12)
plot(sarima_fore, 
     main = "Forecast for the Next 12 Months for Atmospheric Carbon Dioxide", 
     xlab = "Time (Month of Year)", 
     ylab = "Monthly concentration of Atmospheric Carbon Dioxide")
legend("topleft", 
       legend = c("Forecast", "95% Prediction Interval"), 
       col = c("blue", "gray"), 
       lty = c(1, NA), 
       fill = c(NA, "gray"), 
       bty = "n", 
       cex = 0.8) 

par(mfrow = c(1, 1))


par(mfrow = c(1,2))
par(cex.lab = 1, cex.axis = 1, cex.main = 1) 
res_fit3 <- residuals(sarima_model)
Acf(res_fit3, main = "ACF of Residuals(Sarima)", lag.max = 22*12)
pacf(res_fit3, main = "PACF of Residuals(Sarima)",lag.max = 22*12)
par(mfrow = c(1,1))

accuracy_arima_fit1 <- accuracy(sarima_model)


# second arima model (auto arima)
par(mfrow = c(1, 2))
par(cex.lab = 0.7, cex.axis = 0.7, cex.main = 0.7) 

auto.a <- auto.arima(CO2_ts)
auto.a_fit1 <- fitted(auto.a)
summary(auto.a)

par(mfrow = c(1, 2))
par(cex.lab = 0.7, cex.axis = 0.7, cex.main = 0.7) 
plot(CO2_ts, 
     xlab = "Time (Month of Year)", 
     ylab = "Monthly concentration of Atmospheric Carbon Dioxide",
     main = "Time Series with Fitted Values")
lines(auto.a_fit1, col = 2)
legend("topleft", 
       legend = c("Observed", "Fitted (Auto ARIMA)"), 
       col = c("black", "red"), 
       lty = c(1, 1), 
       bty = "n", 
       cex = 0.8)

f2 <- forecast(auto.a, h = 12)
plot(f2, 
     main = "Forecast for the Next 12 Months for Atmospheric Carbon Dioxide", 
     xlab = "Time (Month of Year)", 
     ylab = "Monthly concentration of Atmospheric Carbon Dioxide")
legend("topleft", 
       legend = c("Forecast", "95% Prediction Interval"), 
       col = c("blue", "gray"), 
       lty = c(1, NA), 
       fill = c(NA, "gray"), 
       bty = "n", 
       cex = 0.8)

par(mfrow = c(1, 1))


par(mfrow = c(1,2))
par(cex.lab = 1, cex.axis = 1, cex.main = 1) 
res_fit3 <- residuals(auto.a)
Acf(res_fit3, main = "ACF of Residuals(ARIMA)", lag.max = 22*12)
pacf(res_fit3, main = "PACF of Residuals(ARIMA)",lag.max = 22*12)
par(mfrow = c(1,1))

accuracy_arima_fit2 <- accuracy(f2)


##################################################################################################################
# Methane
##################################################################################################################
df2 <- read.csv("data/global-methane-concentrations.csv")
View(df2)
attach(df2)
head(df2)
Methane = Monthly.concentration.of.atmospheric.methane
Date = Day
tsdisplay(Methane, lag.max = 120)
Methane_ts <- ts(Methane, frequency = 12, start = c(1983, 1))
seasonplot(Methane_ts, ylab = "Monthly concentration of Atmospheric Carbon Dioxide", xlab = "Month of Year",
           main = "Seasonal Plot: concentration of Atmospheric Carbon Dioxide", 
           year.labels = TRUE, year.labels.left = TRUE, 
           col = 1:40, pch = 19)

# Perform STL decomposition
decomposition <- stl(Methane_ts, s.window = 'periodic')

# Extract components into a data frame
components <- data.frame(
  Date = as.Date(time(Methane_ts)), # Convert time series to dates
  Observed = as.numeric(Methane_ts), # Original data
  Seasonal = decomposition$time.series[, "seasonal"],
  Trend = decomposition$time.series[, "trend"],
  Remainder = decomposition$time.series[, "remainder"]
)

# Create individual plots for each component with colors
p1 <- ggplot(components, aes(x = Date, y = Observed)) +
  geom_line(color = "blue") +
  labs(title = "Observed Methane", x = "", y = "Methane") +
  theme_minimal()

p2 <- ggplot(components, aes(x = Date, y = Seasonal)) +
  geom_line(color = "orange") +
  labs(title = "Seasonal Component", x = "", y = "Seasonal") +
  theme_minimal()

p3 <- ggplot(components, aes(x = Date, y = Trend)) +
  geom_line(color = "green") +
  labs(title = "Trend Component", x = "", y = "Trend") +
  theme_minimal()

p4 <- ggplot(components, aes(x = Date, y = Remainder)) +
  geom_line(color = "red") +
  labs(title = "Remainder", x = "Date", y = "Remainder") +
  theme_minimal()

# Arrange plots in a single grid
grid.arrange(p1, p2, p3, p4, ncol = 1)

################################################################################################################
# Linear Regression
################################################################################################################
tt<- 1:NROW(df2)
model1 <- lm(Methane ~ tt)
summary(model1)
anova(model1)
accuracy(model1)

##plot of the model
plot(tt, Methane, xlab="Time(Month of Year)", ylab="Monthly concentration of Atmospheric Carbon Dioxide")
abline(model1, col=6, lwd = 2)

dwtest(model1)

resmodel1<- residuals(model1)
plot(resmodel1,xlab="Time(Month of Year)", ylab="residuals" )
Acf(resmodel1, lag.max = 50*12)
pacf(resmodel1, lag.max = 50*12)
################################################################################
# TSLM model  
################################################################################
ts.plot(Methane_ts, type="o")

## we fit a linear model with the tslm function
TSLM_model2<- tslm(Methane_ts~ trend + season, data = df2)

###obviously it gives the same results of the first model
summary(TSLM_model2)
accuracy(TSLM_model2)

par(mfrow = c(1, 2))
par(cex.lab = 0.7, cex.axis = 0.7, cex.main = 0.6)

plot(Methane_ts,
     xlab = "Time (Month of Year)",
     ylab = "Monthly concentration of Atmospheric Methane",
     main = "Time Series with Fitted Values")
lines(fitted(TSLM_model2), col = 2)


TSLM_fore <- forecast(TSLM_model2, h = 12)
plot(TSLM_fore,
     main = "Forecast for the Next 12 Months for Atmospheric Methane",
     xlab = "Time (Month of Year)",
     ylab = "Monthly concentration of Atmospheric Methane")

par(mfrow = c(1, 1))

TSLM_res<- residuals(TSLM_model2)
par(mfrow = c(1,3))
par(cex.lab = 0.9, cex.axis = 0.9, cex.main = 0.9)
plot(TSLM_res,xlab="Time(Month of Year)", ylab="residuals", main = "TSLM Model Residuals Plot")
Acf(TSLM_res, lag.max = 22*12,xlab="Time(Month of Year)", ylab="ACF for Atmospheric Methane" , main = "ACF for TSLM Model")
pacf(TSLM_res, lag.max = 22*12,xlab="Time(Month of Year)", ylab="PACF for Atmospheric Methane",main = "PACF for TSLM Model")
par(mfrow = c(1,1))
dwtest(TSLM_model2)

################################################################################
##Exponential Smoothing methods
################################################################################
autoplot(Methane_ts)+ylab("Monthly concentration of Atmospheric Methane")+xlab("Time(Month of Year)")

fit1<- ses(Methane_ts, alpha=0.2, initial="simple", h=5)
fit2<- ses(Methane_ts, alpha=0.6, initial="simple", h=5)
fit3<- ses(Methane_ts, h=5)


par(mfrow = c(3,1))

plot(Methane_ts, ylab = "Monthly concentration of Atmospheric Methane",
     xlab = "Time (Month of Year)", main = "Original Time Series with Fit 1")
lines(fitted(fit1), col = "blue", type = "o", lwd = 2)

plot(Methane_ts, ylab = "Monthly concentration of Atmospheric Methane",
     xlab = "Time (Month of Year)", main = "Original Time Series with Fit 2")
lines(fitted(fit2), col = "red", type = "o", lwd = 2)


plot(Methane_ts, ylab = "Monthly concentration of Atmospheric Methane",
     xlab = "Time (Month of Year)", main = "Original Time Series with Fit 3")
lines(fitted(fit3), col = "green", type = "o", lwd = 2)
par(mfrow = c(1,1))


fc<- ses(Methane_ts, h=12)
accuracy(fc)

summary(fc)

autoplot(fc)+
  autolayer(fitted(fc), series="Fitted")+ylab("Monthly concentration of Atmospheric Methane")+xlab("Time (Month of Year)")



##2.Trend methods (Holt method)
fc_Holt<- holt(Methane_ts, h=12)
fc2_Holt<- holt(Methane_ts, damped=T, phi=0.9, h=12)

autoplot(Methane_ts)+
  autolayer(fc_Holt, series="Holt's method", PI=F)+
  autolayer(fc2_Holt, series="Damped Holt's method", PI=F)

###3.Trend and seasonality methods (Holt-Winters method)
autoplot(Methane_ts)

fit1<- hw(Methane_ts, seasonal="additive")
fit2<- hw(Methane_ts, seasonal="multiplicative")

autoplot(Methane_ts, series = "Original Time Series") +
  autolayer(fit1, series = "Additive Holt-Winters", PI = FALSE) +
  autolayer(fit2, series = "Multiplicative Holt-Winters", PI = FALSE) +
  ggtitle("Holt-Winters Forecasts with Different Seasonal Models") +
  xlab("Time (Month of Year)") +
  ylab("Monthly Concentration of Atmospheric Methane") +
  theme_minimal() +
  scale_color_manual(
    values = c("Original Time Series" = "black",
               "Additive Holt-Winters" = "blue",
               "Multiplicative Holt-Winters" = "red")
  ) +
  theme(
    plot.title = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )


# Evaluate performance of fit1 and fit2
accuracy_fit1 <- accuracy(fit1)
accuracy_fit2 <- accuracy(fit2)

par(mfrow = c(2,2))
par(cex.lab = 0.7, cex.axis = 0.7, cex.main = 0.7)
res_fit1 <- residuals(fit1)
Acf(res_fit1, main = "ACF of Residuals (Additive)", lag.max = 22*12)
pacf(res_fit1, main = "PACF of Residuals (Additive)",lag.max = 22*12)

res_fit2 <- residuals(fit2)
Acf(res_fit2, main = "ACF of Residuals (Multiplicative)",lag.max = 22*12)
pacf(res_fit2, main = "PACF of Residuals (Multiplicative)",lag.max = 22*12)
par(mfrow = c(1,1))

########################################################################################################
# SARIMA Models
########################################################################################################
diff1 <- diff(Methane_ts)

# Seasonal differencing to remove the seasonal pattern
diff_seasonal <- diff(diff1, lag = 12)

# Display ACF and PACF after differencing
tsdisplay(diff_seasonal, lag.max = 22*12)


# first arima model
sarima_model <- Arima(Methane_ts, order = c(1, 1, 1), seasonal = list(order = c(0, 1, 0), period = 12))

fit1<- fitted(sarima_model)

par(mfrow = c(1, 2))
par(cex.lab = 0.7, cex.axis = 0.7, cex.main = 0.7)
plot(Methane_ts,
     xlab = "Time (Month of Year)",
     ylab = "Monthly concentration of Atmospheric Methane",
     main = "Time Series with Fitted Values")
lines(fitted(sarima_model), col = 2)
legend("topleft",
       legend = c("Observed", "Fitted (SARIMA)"),
       col = c("black", "red"),
       lty = c(1, 1),
       bty = "n",
       cex = 0.8)


sarima_fore <- forecast(sarima_model, h = 12)
plot(sarima_fore,
     main = "Forecast for the Next 12 Months for Atmospheric Methane",
     xlab = "Time (Month of Year)",
     ylab = "Monthly concentration of Atmospheric Methane")
legend("topleft",
       legend = c("Forecast", "95% Prediction Interval"),
       col = c("blue", "gray"),
       lty = c(1, NA),
       fill = c(NA, "gray"),
       bty = "n",
       cex = 0.8)

par(mfrow = c(1, 1))


par(mfrow = c(1,2))
par(cex.lab = 1, cex.axis = 1, cex.main = 1)
res_fit3 <- residuals(sarima_model)
Acf(res_fit3, main = "ACF of Residuals(Sarima)", lag.max = 22*12)
pacf(res_fit3, main = "PACF of Residuals(Sarima)",lag.max = 22*12)
par(mfrow = c(1,1))

accuracy_arima_fit1 <- accuracy(sarima_model)


# second arima model (auto arima)
par(mfrow = c(1, 2))
par(cex.lab = 0.7, cex.axis = 0.7, cex.main = 0.7)

auto.a <- auto.arima(Methane_ts)
auto.a_fit1 <- fitted(auto.a)
summary(auto.a)

par(mfrow = c(1, 2))
par(cex.lab = 0.7, cex.axis = 0.7, cex.main = 0.7)
plot(Methane_ts,
     xlab = "Time (Month of Year)",
     ylab = "Monthly concentration of Atmospheric Methane",
     main = "Time Series with Fitted Values")
lines(auto.a_fit1, col = 2)
legend("topleft",
       legend = c("Observed", "Fitted (Auto ARIMA)"),
       col = c("black", "red"),
       lty = c(1, 1),
       bty = "n",
       cex = 0.8)

f2 <- forecast(auto.a, h = 12)
plot(f2,
     main = "Forecast for the Next 12 Months for Atmospheric Methane",
     xlab = "Time (Month of Year)",
     ylab = "Monthly concentration of Atmospheric Methane")
legend("topleft",
       legend = c("Forecast", "95% Prediction Interval"),
       col = c("blue", "gray"),
       lty = c(1, NA),
       fill = c(NA, "gray"),
       bty = "n",
       cex = 0.8)

par(mfrow = c(1, 1))


par(mfrow = c(1,2))
par(cex.lab = 1, cex.axis = 1, cex.main = 1)
res_fit3 <- residuals(auto.a)
Acf(res_fit3, main = "ACF of Residuals(ARIMA)", lag.max = 22*12)
pacf(res_fit3, main = "PACF of Residuals(ARIMA)",lag.max = 22*12)
par(mfrow = c(1,1))

accuracy_arima_fit2 <- accuracy(f2)
#########################################################################################################
g1a <- gam(Methane_ts~s(tt))
par(mfrow=c(1,2))
plot(g1a, se=T)
summary(g1a)
