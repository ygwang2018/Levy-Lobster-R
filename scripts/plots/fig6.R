library(dplyr)
library(splines)

test_fig6 <- function() {
  
  # Ensure SEX labels
  if (is.numeric(lobster$SEX)) {
    lobster <- lobster %>%
      mutate(SEX = factor(SEX, levels=c(1,2), labels=c("Female","Male")))
  }
  
  # Split data
  lobster_female <- filter(lobster, SEX == "Female")
  lobster_male   <- filter(lobster, SEX == "Male")
  
  
  # Fit Gamma GLMs
  model_female <- glm(INT ~ ns(PL, df = 3) + ns(NINT, df = 3),
                      data=lobster_female, family=Gamma(link="log"))
  
  model_male <- glm(INT ~ ns(PL, df = 3) + ns(NINT, df = 3),
                    data=lobster_male, family=Gamma(link="log"))
  
  
  # OPEN ONE PNG FOR ALL 8 PANELS
  png("results/figures/fig6.png", width = 1800, height = 2400, res = 150)
  
  # 4 rows × 2 columns = 8 diagnostic plots (Female + Male)
  par(mfrow = c(4, 2))
  
  # FEMALE 
  plot(model_female, which = 1, main = "Residuals vs Fitted - Female")
  plot(model_female, which = 2, main = "Q-Q Plot - Female")
  plot(model_female, which = 3, main = "Scale-Location - Female")
  plot(model_female, which = 4, main = "Residuals vs Leverage - Female")
  
  # MALE
  plot(model_male, which = 1, main = "Residuals vs Fitted - Male")
  plot(model_male, which = 2, main = "Q-Q Plot - Male")
  plot(model_male, which = 3, main = "Scale-Location - Male")
  plot(model_male, which = 4, main = "Residuals vs Leverage - Male")
  
  dev.off()
   
  invisible(list(female = model_female, male = model_male))
}
