install.packages("dplyr")
library(dplyr)

data <- read.csv("C:\\Users\\raulu\\OneDrive\\Documents\\4th Year\\NPS\\Epiphany\\Report\\Datasets\\Master.csv")
data <- na.omit(data)

data_filtered <- data %>%
  select(AMT_INCOME_TOTAL, AMT_CREDIT, EXT_SOURCE_1, EXT_SOURCE_2, EXT_SOURCE_3, DAYS_BIRTH, NAME_INCOME_TYPE)

# Variables: 
# AMT_INCOME_TOTAL – Total annual income of the applicant
# AMT_CREDIT – Total credit amount of the loan.
# EXT_SOURCE_1 – External risk/ credit score from a third-party provider (source 1). - LET'S TAKE A MEAN FROM ALL THE SOURCES!
# DAYS_BIRTH – Age of the client in days (negative values).
# NAME_INCOME_TYPE – Type of income source for the applicant

data_filtered <- data_filtered %>%
  rename(
    Total_Annual_Income = AMT_INCOME_TOTAL,
    Total_Loan_Credit = AMT_CREDIT,
    Income_Source = NAME_INCOME_TYPE
  ) %>%
  mutate(
    External_Credit_Score = rowMeans(select(., EXT_SOURCE_1, EXT_SOURCE_2, EXT_SOURCE_3), na.rm = TRUE),
    Age = abs(DAYS_BIRTH) / 365
  ) %>%
  select(-DAYS_BIRTH, -EXT_SOURCE_1, -EXT_SOURCE_2, -EXT_SOURCE_3)  # Remove old columns

head(data_filtered)

# Reset row names
rownames(data_filtered) <- NULL

# View the first few rows again
head(data_filtered)
# Define the target median value for Age
target_age <- 39.48

# Sort the dataset by Age
data <- data_filtered %>%
  arrange(Age)

# Get the 50 closest values below 39.48
below_age <- data %>%
  filter(Age < target_age) %>%
  slice_tail(n = 50)

# Get the 50 closest values above 39.48
above_age <- data %>%
  filter(Age >= target_age) %>%
  slice_head(n = 50)

# Combine both sets
data <- bind_rows(below_age, above_age)

# View the first few rows to confirm
head(data)

# Save the filtered dataset
write.csv(data, "C:\\Users\\raulu\\OneDrive\\Documents\\4th Year\\NPS\\Epiphany\\Report\\final_data.csv", row.names = FALSE)

hist(data$Age, freq = FALSE)
hist(data$External_Credit_Score, freq = FALSE)
# hist(data_filtered_age$Total_Annual_Income)

install.packages("KernSmooth")
install.packages("locpol")
require(KernSmooth)
require(locpol)

data <- data_filtered

dpill(data$Age, data$External_Credit_Score) # 0.05557425
Age.grid <- seq(20, 70, by=0.5)
fit<- locpol(External_Credit_Score~Age, bw=3.48996, kernel=gaussK, xeval=Age.grid, data=data)
plot(fit)

# Yes this is now a bit better. Note that the difference between the optimal bandwidths of the two approaches is larger for the
# lidar data, which are heteroscedastic, and hence may more likely deliver a poorly performing ROT bandwidth.