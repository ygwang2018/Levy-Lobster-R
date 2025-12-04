
library(readr)
library(dplyr)
library(splines)

# If SEX is numeric 1/2, first make it labelled:
# lobster <- lobster |>
#   mutate(SEX = factor(SEX, levels = c(1, 2), labels = c("Female", "Male")))

# Split by sex using factor labels
lobster_female <- filter(lobster, SEX == "Female")
lobster_male   <- filter(lobster, SEX == "Male")

# Fit Gamma GLMs with spline terms
model_female <- glm(
  INT ~ ns(PL, df = 3) + ns(NINT, df = 3),
  data   = lobster_female,
  family = Gamma(link = "log")
)

model_male <- glm(
  INT ~ ns(PL, df = 3) + ns(NINT, df = 3),
  data   = lobster_male,
  family = Gamma(link = "log")
)

# Open a PNG device
png("results/figures/diagnostic_plots_female_male.png", 
    width = 1200, height = 800, res = 150)

# --- FEMALE ----
par(mfrow = c(2, 2))
plot(model_female, which = 1,
     main = "Residuals vs Fitted - Female")
plot(model_female, which = 2,
     main = "QQ Plot - Female")
plot(model_female, which = 3,
     main = "Scale-Location - Female")
plot(model_female, which = 4,
     main = "Residuals vs Leverage - Female")

# --- MALE ----
par(mfrow = c(2, 2))
plot(model_male, which = 1,
     main = "Residuals vs Fitted - Male")
plot(model_male, which = 2,
     main = "QQ Plot - Male")
plot(model_male, which = 3,
     main = "Scale-Location - Male")
plot(model_male, which = 4,
     main = "Residuals vs Leverage - Male")

# Close the device (important!)
dev.off()
