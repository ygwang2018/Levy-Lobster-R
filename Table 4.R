library(dplyr)

# Prepare data
lobster$SEX <- factor(lobster$SEX, levels=c(1,2), labels=c("Female","Male"))

# Function to fit Gamma GLM and extract formatted table row
fit_ip_model <- function(data, sex_label) {

  # Fit Gamma GLM for INT model
  fit_obj <- glm(
    INT ~ PL + NINT,
    data = data,
    family = Gamma(link="log")
  )

  # Extract coefficients and SE
  coefs <- summary(fit_obj)$coefficients

  b0    <- coefs["(Intercept)", 1]
  b0_se <- coefs["(Intercept)", 2]

  b1    <- coefs["PL", 1]
  b1_se <- coefs["PL", 2]

  b2    <- coefs["NINT", 1]
  b2_se <- coefs["NINT", 2]

  # Dispersion parameter φ
  phi <- summary(fit_obj)$dispersion

  # Return formatted output
  data.frame(
    Sex = sex_label,
    Intercept = sprintf("%.3f (%.3f)", b0, b0_se),
    Premoult_length = sprintf("%.3f (%.3f)", b1, b1_se),
    Previous_IP = sprintf("%.3f (%.3f)", b2, b2_se),
    Dispersion = sprintf("%.3f", phi),
    stringsAsFactors = FALSE
  )
}

# Fit models by sex
female_table <- fit_ip_model(filter(lobster, SEX=="Female"), "Female")
male_table   <- fit_ip_model(filter(lobster, SEX=="Male"),   "Male")

# Combine final table
final_table <- rbind(female_table, male_table)

final_table
