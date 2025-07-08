# Load your dataset
lobster <- read_csv("data/lobster.csv")
# Clean the data
lobster <- lobster %>%
  filter(!is.na(INT), INT > 0, !is.na(PL), !is.na(NINT), !is.na(SEX))
# Function to fit model and plot diagnostics
plot_diagnostics_by_sex <- function(data, sex_label) {
  # Fit Gamma GLM with spline terms
  model <- glm(INT ~ ns(PL, df = 3) + ns(NINT, df = 3),
               data = data, family = Gamma(link = "log"))
  # Set up plotting area
  par(mfrow = c(2, 2))
  plot(model, which = 1, main = paste("Residuals vs Fitted -", sex_label))
  plot(model, which = 2, main = paste("QQ Plot -", sex_label))
  plot(model, which = 3, main = paste("Scale-Location -", sex_label))
  plot(model, which = 4, main = paste("Residuals vs Leverage -", sex_label))
  # Optionally return model if needed
  return(model)
}
# Fit and plot for females
lobster_female <- filter(lobster, SEX == 1)
model_female <- plot_diagnostics_by_sex(lobster_female, "Female")
# Fit and plot for males
lobster_male <- filter(lobster, SEX == 2)
model_male <- plot_diagnostics_by_sex(lobster_male, "Male")