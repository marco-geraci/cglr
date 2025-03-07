rm(list = ls())

source("cglr.R")
sourceCpp("cglr.cpp")

# Read data from csv 
# Data available at https://www.kaggle.com/datasets/berkerisen/wind-turbine-scada-dataset

dd <- read.csv("turbine.csv")

# Convert the string to a date-time format
dd$date_time <- as.POSIXct(dd$date, format = "%d %m %Y %H:%M")

# Remove missing (DST?)
dd <- subset(dd, !is.na(dd$date_time))

# Extract the date
dd$date <- format(dd$date_time, "%Y-%m-%d")

# Extract the month
dd$month <- factor(as.numeric(format(dd$date_time, "%m")), levels = 1:12, labels = month.abb)

# Extract hours
dd$hour <- as.numeric(format(dd$date_time, "%H"))

# Convert degrees to radians
dd$theta <- dd$angle/180*pi

# Take a sample
n <- nrow(dd)
set.seed(123)
dds <- dd[sample(1:n, 5000),]

#########################################################
#########################################################

Hv <- 20 # number of quadrature points

# Fit CGLR model with 1 linear outcome
fit0_cglr <- cglr(speed|theta ~ 1, data = dds, control = cglrControl(optimizer = "optim", Hv = Hv, verbose = TRUE))
fit1_cglr <- cglr(speed|theta ~ month, data = dds, control = cglrControl(optimizer = "optim", Hv = Hv))
fit2_cglr <- cglr(speed|theta ~ month, mu.formula = ~ month, alpha.formula = ~ 1, data = dds, control = cglrControl(optimizer = "optim", Hv = Hv))

# Fit CPNR model
fit0_cpnr <- cpnr(speed|theta ~ 1, data = dds, control = cglrControl(optimizer = "optim", Hv = Hv))
fit1_cpnr <- cpnr(speed|theta ~ month, data = dds, control = cglrControl(optimizer = "optim", Hv = Hv))

#########################################################
#########################################################

# Plot fit0
library(ggplot2)

nk <- 200
tt <- seq(0, 2*pi, length = nk)
ww <- seq(0, 25, length = nk)
df <- data.frame(theta = tt, speed = ww)

df$value_cglr <- dtheta.cglr(fit0_cglr, data = df)
df$value_cpnr <- dtheta.cpnr(fit0_cpnr, data = df)

bks <- c(0, pi/2, pi, 3/2*pi, 2*pi)
xlabs <- bks*180/pi

ggplot(dds, aes(x = theta)) + geom_histogram(aes(y = after_stat(density)), alpha=0.6, bins = 30) + geom_line(data = df, aes(x = theta, y = value_cglr), linewidth = 1.2)  + geom_line(data = df, aes(x = theta, y = value_cpnr), linewidth = 1.2, linetype = "dashed") + scale_x_continuous("Wind direction", breaks = bks, labels = xlabs) + theme(legend.position="none", text = element_text(size=10)) + ylab("Density")

df$value_cglr <- dw.cglr(fit0_cglr, variable = "speed", data = df)
df$value_cpnr <- dw.cpnr(fit0_cpnr, variable = "speed", data = df)

bks <- seq(0, 25, by = 5)
xlabs <- bks

ggplot(dds, aes(x = speed)) + geom_histogram(aes(y = after_stat(density)), alpha=0.6, bins = 30) + geom_line(data = df, aes(x = speed, y = value_cglr), linewidth = 1.2) + geom_line(data = df, aes(x = speed, y = value_cpnr), linewidth = 1.2, linetype = "dashed") + scale_x_continuous("Wind speed (m/s)", breaks = bks, labels = xlabs) + theme(legend.position="none", text = element_text(size=10)) + ylab("Density")

# Plot fit1

nk <- 200
tt <- seq(0, 2*pi, length = nk)
ww <- seq(0, 25, length = nk)
df <- expand.grid(theta = tt, month = levels(dd$month))
df$speed <- ww
df$value_cglr <- dtheta.cglr(fit1_cglr, data = df)
df$value_cpnr <- dtheta.cpnr(fit1_cpnr, data = df)

ggplot(dds, aes(x = theta, group = month)) + geom_histogram(aes(y = after_stat(density)), alpha=0.6, bins = 30) + facet_wrap( ~ month, ncol = 3, scales = "free_y") + geom_line(data = df, aes(x = theta, y = value_cglr), linewidth = 1.2) + geom_line(data = df, aes(x = theta, y = value_cpnr), linewidth = 1.2, linetype = "dashed") + scale_x_continuous("Wind direction", breaks = bks, labels = xlabs) + theme(legend.position="none", text = element_text(size=10)) + ylab("Density")

df$value_cglr <- dw.cglr(fit1_cglr, variable = "speed", data = df)
df$value_cpnr <- dw.cpnr(fit1_cpnr, variable = "speed", data = df)

bks <- seq(0, 25, by = 5)
xlabs <- bks

ggplot(dds, aes(x = speed, group = month)) + geom_histogram(aes(y = after_stat(density)), alpha=0.6, bins = 30) + facet_wrap( ~ month, ncol = 3, scales = "fixed") + geom_line(data = df, aes(x = speed, y = value_cglr), linewidth = 1.2) + geom_line(data = df, aes(x = speed, y = value_cpnr), linewidth = 1.2, linetype = "dashed") + scale_x_continuous("Wind speed (m/s)", breaks = bks, labels = xlabs) + theme(legend.position="none", text = element_text(size=10)) + ylab("Density")


# Plot fit2

nk <- 200
tt <- seq(0, 2*pi, length = nk)
ww <- seq(0, 25, length = nk)
df <- expand.grid(theta = tt, month = levels(dd$month))
df$speed <- ww
df$value_cglr <- dtheta.cglr(fit2_cglr, data = df)
df$value_cpnr <- dtheta.cpnr(fit1_cpnr, data = df)

bks <- c(0, pi/2, pi, 3/2*pi, 2*pi)
xlabs <- bks*180/pi

ggplot(dds, aes(x = theta, group = month)) + geom_histogram(aes(y = after_stat(density)), alpha=0.6, bins = 30) + facet_wrap( ~ month, ncol = 3, scales = "free_y") + geom_line(data = df, aes(x = theta, y = value_cglr), linewidth = 1.2) + geom_line(data = df, aes(x = theta, y = value_cpnr), linewidth = 1.2, linetype = "dashed") + scale_x_continuous("Wind direction", breaks = bks, labels = xlabs) + theme(legend.position="none", text = element_text(size=10)) + ylab("Density")

df$value_cglr <- dw.cglr(fit2_cglr, variable = "speed", data = df)
df$value_cpnr <- dw.cpnr(fit1_cpnr, variable = "speed", data = df)

bks <- seq(0, 25, by = 5)
xlabs <- bks

ggplot(dds, aes(x = speed, group = month)) + geom_histogram(aes(y = after_stat(density)), alpha=0.6, bins = 30) + facet_wrap( ~ month, ncol = 3, scales = "fixed") + geom_line(data = df, aes(x = speed, y = value_cglr), linewidth = 1.2) + geom_line(data = df, aes(x = speed, y = value_cpnr), linewidth = 1.2, linetype = "dashed") + scale_x_continuous("Wind speed (m/s)", breaks = bks, labels = xlabs) + theme(legend.position="none", text = element_text(size=10)) + ylab("Density")

