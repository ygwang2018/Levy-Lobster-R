library(dplyr)
library(splines)

test_fig7 <- function() {
  cat("Running test_fig7...\n")
  
  # Ensure SEX has correct labels
  if (is.numeric(lobster$SEX)) {
    lobster <- lobster %>%
      mutate(SEX = factor(SEX, levels=c(1,2), labels=c("Female","Male")))
  }
  
  # Split data
  lobster_female <- filter(lobster, SEX == "Female")
  lobster_male   <- filter(lobster, SEX == "Male")
  
  if (nrow(lobster_female) == 0 | nrow(lobster_male) == 0) {
    stop("ERROR: No rows found for one or both sexes. Check SEX coding in `lobster`.")
  }
  
  # Fit Gamma GLMs
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
  
  # Create output folder
  if (!dir.exists("results/figures")) dir.create("results/figures", recursive = TRUE)
  
  # Open PNG device (big canvas)
  png("results/figures/fig7.png", width = 1800, height = 1200, res = 150)
  
  # ---------------------------
  # PAGE 1 — FEMALE
  # ---------------------------
  par(mfrow = c(2, 2))
  plot(model_female, which = 1, main = "Residuals vs Fitted - Female")
  plot(model_female, which = 2, main = "Q-Q Plot - Female")
  plot(model_female, which = 3, main = "Scale-Location - Female")
  plot(model_female, which = 4, main = "Residuals vs Leverage - Female")
  
  # Force new page
  dev.flush()
  
  # ---------------------------
  # PAGE 2 — MALE
  # ---------------------------
  par(mfrow = c(2, 2))
  plot(model_male, which = 1, main = "Residuals vs Fitted - Male")
  plot(model_male, which = 2, main = "Q-Q Plot - Male")
  plot(model_male, which = 3, main = "Scale-Location - Male")
  plot(model_male, which = 4, main = "Residuals vs Leverage - Male")
  
  dev.off()
  
  cat("fig7 saved to results/figures/fig7.png\n")
  invisible(list(female = model_female, male = model_male))
}
