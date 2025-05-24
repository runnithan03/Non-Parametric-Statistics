install.packages(c("quantmod", "boot", "bayesboot", "BNPmix", "ggplot2", "reshape2", "dplyr", "xgboost", "Matrix", "data.table", "caret", "MLmetrics","coda","MCMCpack","moments"))

library(quantmod)
library(boot)
library(bayesboot)
library(BNPmix)
library(ggplot2)
library(reshape2)
library(dplyr)
library(xgboost)
library(Matrix)
library(data.table)
library(caret)
library(MLmetrics)
library(coda)
library(MCMCpack)
library(moments)

vix_data <- read.csv("C:\\Users\\raulu\\OneDrive\\Documents\\4th Year\\NPS\\Michaelmas\\Report\\Raul_Unnithan_Data.csv")

# Dataset Pre-Processing
vix_data_close <- vix_data[, c("Date", "Close")]
vix_data_close <- na.omit(vix_data_close)
vix_data_close$Date <- as.Date(vix_data_close$Date)
vix_data_close$Close <- as.numeric(vix_data_close$Close)

# Initial Time Series
ggplot(vix_data_close, aes(x = Date, y = Close)) +
  geom_line(color = "blue") +
  labs(title = "Time-Series of VIX Closing Prices",
       x = "Date",
       y = "Closing Price") +
  theme_minimal()

# Q-Q plot for close price
close_price <- vix_data_close$Close  
qqnorm(close_price, main = "QQ Plot of Close Price")
qqline(close_price, col = "blue")

# Compute Log Returns (Percentage Change)
vix_returns <- diff(log(vix_data_close$Close)) * 100
cat("Summary of VIX Returns:\n")
print(summary(vix_returns))

# Classical Bootstrap Analysis
stat_mean <- function(data, indices) {
  mean(data[indices], na.rm = TRUE)
}

stat_variance <- function(data, indices) {
  var(data[indices], na.rm = TRUE)
}

set.seed(123)
boot_mean <- boot(data = vix_returns, statistic = stat_mean, R = 2000)
boot_variance <- boot(data = vix_returns, statistic = stat_variance, R = 2000)

ci_mean <- boot.ci(boot_mean, type = "bca")
ci_variance <- boot.ci(boot_variance, type = "bca")

cat("Classical Bootstrap Results:\n")
print(ci_mean)
print(ci_variance)

# Plot Bootstrap Distributions
ggplot(data.frame(boot_mean = boot_mean$t), aes(x = boot_mean)) +
  geom_histogram(fill = "green", alpha = 0.5, bins = 30, color = "black") +
  labs(title = "Bootstrap Distribution for Mean",
       x = "Bootstrap Mean", y = "Frequency") +
  theme_minimal()

ggplot(data.frame(boot_variance = boot_variance$t), aes(x = boot_variance)) +
  geom_histogram(fill = "purple", alpha = 0.5, bins = 30, color = "black") +
  labs(title = "Bootstrap Distribution for Variance",
       x = "Bootstrap Variance", y = "Frequency") +
  theme_minimal()

# Bayesian Bootstrap Analysis
set.seed(123)
bayes_boot_median <- bayesboot(vix_returns, statistic = median, R = 2000)
bayes_boot_mean <- bayesboot(vix_returns, statistic = weighted.mean, R = 2000, use.weights = TRUE)

ggplot() +
  geom_density(aes(bayes_boot_median$V1), fill = "blue", alpha = 0.5) +
  geom_density(aes(bayes_boot_mean$V1), fill = "red", alpha = 0.5) +
  labs(title = "Bayesian Bootstrap Posterior Distributions",
       x = "Statistic", y = "Density") +
  theme_minimal()

# Compare Classical and Bayesian Bootstrap (Overlay Plot)
bootstrap_data <- data.frame(
  Classical = scale(boot_mean$t),
  Bayesian = scale(bayes_boot_mean$V1)
)

# Convert bootstrap_data from a wide to long format
bootstrap_melted <- melt(bootstrap_data) 

ggplot(bootstrap_melted, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(
    values = c("Classical" = "green", "Bayesian" = "blue"),
    name = "Bootstrap Method:"
  ) +
  labs(
    title = "Overlay of Classical vs Bayesian Bootstrap for Mean (Normalised)",
    x = "Normalised Bootstrap Statistic",
    y = "Density"
  ) +
  theme_minimal()

# Bayesian Nonparametric Analysis (using DPMM)
set.seed(123)
dpmm_result <- PYdensity(
  y = vix_returns,
  mcmc = list(niter = 2000, nburn = 500),
  prior = list(strength = 1, discount = 0, hyper = FALSE),
  output = list(grid = seq(min(vix_returns), max(vix_returns), length.out = 100), out_param = TRUE)
)

cat("DPMM Summary:\n")
print(summary(dpmm_result))

plot(dpmm_result, band = TRUE, show_hist = TRUE, 
     xlab = "VIX Returns", ylab = "Density", main = "DPMM Posterior Density")

# Trace plot of cluster counts over iterations
cluster_counts <- apply(dpmm_result$clust, 1, function(x) length(unique(x))) # Summarise cluster counts over iterations
cluster_counts_mcmc <- as.mcmc(cluster_counts) # Convert to mcmc object for trace plot
traceplot(cluster_counts_mcmc, main = "Trace Plot of Cluster Counts Over Iterations", ylab = "Number of Clusters", xlab = "Iterations")

# Cluster Analysis based on MCMC iterations
cluster_assignment <- apply(dpmm_result$clust, 1, function(x) {
  valid_clusters <- x[x > 0]
  if (length(valid_clusters) > 0) {
    as.numeric(names(sort(table(valid_clusters), decreasing = TRUE))[1])
  } else {
    NA
  }
})

table(cluster_assignment)

vix_clustered <- data.frame(
  Date = tail(vix_data_close$Date, length(cluster_assignment)),
  Returns = tail(vix_returns, length(cluster_assignment)),
  Cluster = factor(cluster_assignment)
) %>%
  na.omit() %>%
  filter(Cluster %in% names(table(Cluster)[table(Cluster) >= 20])) # keeps only clusters with >= 20 data points

ggplot(vix_clustered, aes(x = Date, y = Returns, color = Cluster)) +
  geom_line() +
  facet_wrap(~ Cluster, scales = "free_y") +
  theme_minimal()

# Mean density across MCMC iterations
mean_density <- colMeans(dpmm_result$density)

# Grid value corresponding to the peak density
max_density_index <- which.max(mean_density)
peak_value <- dpmm_result$grideval[max_density_index]

cat("Peak occurs at VIX Returns:", peak_value, "\n")

# Analysing the peak to find the most contributing regimes
vix_clustered$Date <- as.Date(vix_clustered$Date, origin = "1970-01-01")

peak_regime <- vix_clustered %>%
  filter(Returns >= (peak_value - 0.2) & Returns <= (peak_value + 0.2)) %>%
  count(Cluster, sort = TRUE)

cat("Regime contributing to the peak:\n")
print(peak_regime) # cluster 3 has the most observations

# Cluster 3 analysis
cluster_3_stats <- vix_clustered %>%
  filter(Cluster == 3) %>%
  summarise(
    Mean = mean(Returns),
    Variance = var(Returns),
    Skewness = skewness(Returns),
    Count = n()
  )

# XGBoost Model for Return Prediction (with RSR Calculation)
vix_clustered <- vix_clustered %>%
  mutate(Cluster = as.numeric(Cluster))

set.seed(123)
train_index <- sample(1:nrow(vix_clustered), size = 0.8 * nrow(vix_clustered))
train_data <- vix_clustered[train_index, ]
test_data <- vix_clustered[-train_index, ]

dtrain <- xgb.DMatrix(data = as.matrix(train_data[, c("Returns", "Cluster")]), label = train_data$Returns)
dtest <- xgb.DMatrix(data = as.matrix(test_data[, c("Returns", "Cluster")]))

xgb_model <- xgb.train(
  params = list(objective = "reg:squarederror", 
                eta = 0.01, 
                max_depth = 6,
                subsample = 0.8, 
                colsample_bytree = 0.8),
  data = dtrain,
  nrounds = 500
)

# Predictions and MSE
predictions <- predict(xgb_model, dtest)
mse <- mean((predictions - test_data$Returns)^2)

# RMSE and RSR Calculation
rmse <- sqrt(mse)
std_actual <- sd(test_data$Returns)
rsr <- rmse / std_actual

cat("MSE:", mse, "\n")
cat("RMSE:", rmse, "\n")
cat("RSR (Relative Standardised Residual):", rsr, "\n")

vix_all<- vix_data[, c("Date", "Open", "Close")]
vix_all <- na.omit(vix_all)
vix_all$Date <- as.Date(vix_all$Date)
vix_all$Open <- as.numeric(vix_all$Open)
vix_all$Close <- as.numeric(vix_all$Close)

# Compute Log Returns (as percentage changes) for Open and Close prices
log_returns_open <- diff(log(vix_all$Open)) * 100  
log_returns_close <- diff(log(vix_all$Close)) * 100

# Save the log returns as a new dataframe
vix_returns_all <- data.frame(
  Date = vix_all$Date[-1],  # Exclude the first date since there's no log return for it
  Log_Open = log_returns_open,
  Log_Close = log_returns_close
)

# Display the resulting dataframe
print(vix_returns_all)

cat("Summary of VIX Returns:\n")
print(summary(vix_returns_all))

# Response and predictor selection
Log_Close <- vix_returns_all$Log_Close
Log_Open <- vix_returns_all$Log_Open

# Calibration of the Pitman-Yor prior
n <- length(Log_Open)
prior <- PYcalibrate(Ek = 3, n = n, discount = 0)
print(prior)

# Define MCMC settings
mcmc <- list(niter = 11000, nburn = 1000)

# Define output grid for log_open and log_close
grid_y <- seq(min(Log_Open), max(Log_Open), length.out = 100)
grid_x <- quantile(Log_Close, probs = c(0.01, 0.25, 0.50, 0.75, 0.99))

set.seed(123)
# Fit the PYregression model
fit_reg <- PYregression(
  y = Log_Open,
  x = Log_Close,
  prior = list(strength = prior$strength, discount = 0, hyper = FALSE),
  mcmc = mcmc,
  output = list(grid_y = grid_y, grid_x = grid_x, out_type = "FULL", out_param = TRUE)
)

summary(fit_reg)

# Summarise clusters
cluster_labels <- fit_reg$clust
num_clusters <- apply(cluster_labels, 1, function(row) length(unique(row)))
cluster_counts <- table(num_clusters)
cluster_counts/sum(cluster_counts)

# Cluster probabilities Plot
barplot(cluster_counts / sum(cluster_counts), col = "lightblue",
        main = "Cluster Distribution for Open Return", xlab = "Number of Clusters", ylab = "Probability")

# Prepare regression plot data
regplot <- data.frame(
  dens = as.vector(apply(fit_reg$density, c(1, 2), mean)),
  qlow = as.vector(apply(fit_reg$density, c(1, 2), quantile, probs = 0.025)),
  qupp = as.vector(apply(fit_reg$density, c(1, 2), quantile, probs = 0.975)),
  grid = rep(grid_y, length(grid_x)),
  label = factor(rep(paste("Open =", round(grid_x, 1)), each = length(grid_y)),
                 levels = paste("Open =", round(grid_x, 1)))
)

# Conditional Density plot
ggplot(regplot) +
  theme_bw() +
  geom_line(aes(x = grid, y = dens), color = "blue") +
  geom_ribbon(aes(x = grid, ymin = qlow, ymax = qupp), fill = "lightblue", alpha = 0.3) +
  facet_wrap(~label, ncol = 3, scales = "free_y") +
  labs(x = "Close Return", y = "Density", title = "Posterior Conditional Densities of Close Return Given Open Return")

