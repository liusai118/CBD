library(forecast)
hyper <- read.csv("IHME-GBD_2021_DATA-be1250ce-1.csv")
hyper <- hyper[hyper$metric == "1",]
####
mental <- read.csv("IHME-GBD_2021_DATA-066dd0fd-1.csv")
mental_1 <- mental [mental$metric == "1" & mental$measure=="1",]
mental_2 <- mental [mental$metric == "1" & mental$measure=="2",]
mental_5 <- mental [mental$metric == "1" & mental$measure=="5",]
######

mental <- read.csv("IHME-GBD_2021_DATA-9bfe0a1f-1.csv")
mental_1 <- mental [mental$metric == "1" & mental$measure=="1",]
mental_2 <- mental [mental$metric == "1" & mental$measure=="2",]
mental_5 <- mental [mental$metric == "1" & mental$measure=="5",]
########
mental <- read.csv("IHME-GBD_2021_DATA-be1250ce-1.csv")
mental_1 <- mental [mental$metric == "1" & mental$measure=="1",]
mental_2 <- mental [mental$metric == "1" & mental$measure=="2",]
mental_5 <- mental [mental$metric == "1" & mental$measure=="5",]

########
mental <- read.csv("IHME-GBD_2021_DATA-dd89b590-1.csv")
mental_1 <- mental [mental$metric == "1" & mental$measure=="1",]
mental_2 <- mental [mental$metric == "1" & mental$measure=="2",]
mental_5 <- mental [mental$metric == "1" & mental$measure=="5",]
mental_5 <- mental_5[order(mental_5$year),]
human <- read.csv("IHME-GBD_2021_DATA-8f4dc684-1.csv")
ts_data <- ts(mental_5$val, start = min(mental_5$year), frequency = 2)


model <- auto.arima(ts_data)


future_forecast <- forecast(model, h = 28, level = 95)


future_years <- seq(max(mental_5$year) + 1, max(mental_5$year) + 28, by = 1)
future_data <- data.frame(
  year = future_years,
  val = future_forecast$mean,
  lower = future_forecast$lower[, 1], 
  upper = future_forecast$upper[, 1]
)

last_point <- mental_5[nrow(mental_5), c("year", "val", "lower", "upper")]
future_data <- rbind(last_point, future_data)


full_data <- rbind(
  mental_5[, c("year", "val", "lower", "upper")],
  future_data
)
full_data<- full_data[!duplicated(full_data$year),]

ggplot(full_data, aes(x = year)) +

  geom_ribbon(data = mental_5, aes(ymin = lower, ymax = upper), fill = "#B0B0DA", alpha = 0.3) +
  geom_line(data = mental_5, aes(y = val), color = "#B0B0DA", size = 1) +
  # geom_ribbon(data = future_data, aes(ymin = lower, ymax = upper), fill = "#D2D2EA", alpha = 0.3) +
  geom_line(data = future_data, aes(y = val), color = "#D2D2EA", size = 1) +
  # 
  labs(title = "区间折线图") +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )
##################
library(forecast)
library(ggplot2)


mental_5 <- human
mental_5 <- mental_5[order(mental_5$year),]
ts_data <- ts(mental_5$val, start = min(mental_5$year), frequency = 2)


model_ets <- ets(ts_data)
future_forecast <- forecast(model_ets, h = 28, level = 90)


future_years <- seq(max(mental_5$year) + 1, max(mental_5$year) + 28, by = 1)
future_data <- data.frame(
  year = future_years,
  val = future_forecast$mean,
  lower = future_forecast$lower[, 1],
  upper = future_forecast$upper[, 1]
)


last_point <- mental_5[nrow(mental_5), c("year", "val", "lower", "upper")]
future_data <- rbind(last_point, future_data)


full_data <- rbind(
  mental_5[, c("year", "val", "lower", "upper")],
  future_data
)
full_data <- full_data[!duplicated(full_data$year),]

# 绘图
ggplot(full_data, aes(x = year)) +

  geom_ribbon(data = mental_5, aes(ymin = lower, ymax = upper), fill = "#E1A09B", alpha = 0.3) +
  geom_line(data = mental_5, aes(y = val), color = "#E1A09B", size = 1) +
  geom_ribbon(data = future_data, aes(ymin = lower, ymax = upper), fill = "#EABCB8", alpha = 0.3) +
  geom_line(data = future_data, aes(y = val), color = "#EABCB8", size = 1) +
  labs(title = "区间折线图") +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )
############
mental_5 <- mental_5[order(mental_5$year),]
ts_data <- ts(mental_5$val, start = min(mental_5$year), frequency = 2)
model <- auto.arima(ts_data)
future_forecast <- forecast(model, h = 28, level = 99)

future_years <- seq(max(mental_5$year) + 1, max(mental_5$year) + 28, by = 1)
future_data <- data.frame(
  year = future_years,
  val = future_forecast$mean,
  lower = future_forecast$lower[, 1],  
  upper = future_forecast$upper[, 1]
)


last_point <- mental_5[nrow(mental_5), c("year", "val", "lower", "upper")]
future_data <- rbind(last_point, future_data)


full_data <- rbind(
  mental_5[, c("year", "val", "lower", "upper")],
  future_data
)
full_data<- full_data[!duplicated(full_data$year),]
# 绘图
ggplot(full_data, aes(x = year)) +

  geom_ribbon(data = mental_5, aes(ymin = lower, ymax = upper), fill = "#7DAFD1", alpha = 0.6) +
  geom_line(data = mental_5, aes(y = val), color = "#7DAFD1", size = 1) +

  geom_ribbon(data = future_data, aes(ymin = lower, ymax = upper), fill = "#D1E3EF", alpha = 0.3) +
  geom_line(data = future_data, aes(y = val), color = "#D1E3EF", size = 1) +

  labs(title = "区间折线图") +

  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )
###########
mental_5 <- mental_5[order(mental_5$year),]
ts_data <- ts(mental_5$val, start = min(mental_5$year), frequency = 2)


model <- auto.arima(ts_data)

future_forecast <- forecast(model, h = 28, level = 99)


future_years <- seq(max(mental_5$year) + 1, max(mental_5$year) + 28, by = 1)
future_data <- data.frame(
  year = future_years,
  val = future_forecast$mean,
  lower = future_forecast$lower[, 1],  
  upper = future_forecast$upper[, 1]
)


last_point <- mental_5[nrow(mental_5), c("year", "val", "lower", "upper")]
future_data <- rbind(last_point, future_data)

full_data <- rbind(
  mental_5[, c("year", "val", "lower", "upper")],
  future_data
)
full_data<- full_data[!duplicated(full_data$year),]

ggplot(full_data, aes(x = year)) +

  geom_ribbon(data = mental_5, aes(ymin = lower, ymax = upper), fill = "#C7690B", alpha = 0.6) +
  geom_line(data = mental_5, aes(y = val), color = "#C7690B", size = 1) +

  geom_ribbon(data = future_data, aes(ymin = lower, ymax = upper), fill = "#F8BE84", alpha = 0.3) +
  geom_line(data = future_data, aes(y = val), color = "#F8BE84", size = 1) +

  labs(title = "区间折线图") +

  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )
