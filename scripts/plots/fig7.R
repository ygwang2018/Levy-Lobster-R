library(readr)
library(dplyr)
library(splines)

test_fig7 <- function() {
  cat("Running test for fig7...\n")
  
  # Split by sex (assumes lobster already loaded and SEX is factor with labels)
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
  if (!dir.exists("results/figures")) dir.create("results/figures", recursive = TRUE)
  png("results/figures/fig7.png", width = 1200, height = 800, res = 150)
  
  # --- FEMALE ----
  par(mfrow = c(2, 2))
  plot(model_female, which = 1, main = "Residuals vs Fitted - Female")
  plot(model_female, which = 2, main = "QQ Plot - Female")
  plot(model_female, which = 3, main = "Scale-Location - Female")
  plot(model_female, which = 4, main = "Residuals vs Leverage - Female")
  
  # --- MALE ----
  par(mfrow = c(2, 2))
  plot(model_male, which = 1, main = "Residuals vs Fitted - Male")
  plot(model_male, which = 2, main = "QQ Plot - Male")
  plot(model_male, which = 3, main = "Scale-Location - Male")
  plot(model_male, which = 4, main = "Residuals vs Leverage - Male")
  
  # Close the device
  dev.off()
  
  invisible(list(female = model_female, male = model_male))
}
