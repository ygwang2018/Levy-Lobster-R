test_table4 <- function() {
  library(dplyr)
  
  # Ensure SEX is labeled
  lobster$SEX <- factor(lobster$SEX, levels = c(1,2), labels = c("Female","Male"))
  
  # Helper to fit Gamma GLM and extract formatted row
  fit_ip_model <- function(data, sex_label) {
    fit_obj <- glm(
      INT ~ PL + NINT,
      data = data,
      family = Gamma(link = "log")
    )
    coefs <- summary(fit_obj)$coefficients
    phi   <- summary(fit_obj)$dispersion
    
    data.frame(
      Sex = sex_label,
      Intercept        = sprintf("%.3f (%.3f)", coefs["(Intercept)",1], coefs["(Intercept)",2]),
      Premoult_length  = sprintf("%.3f (%.3f)", coefs["PL",1], coefs["PL",2]),
      Previous_IP      = sprintf("%.3f (%.3f)", coefs["NINT",1], coefs["NINT",2]),
      Dispersion       = sprintf("%.3f", phi),
      stringsAsFactors = FALSE
    )
  }
  
  # Fit models by sex
  female_table <- fit_ip_model(filter(lobster, SEX == "Female"), "Female")
  male_table   <- fit_ip_model(filter(lobster, SEX == "Male"),   "Male")
  
  # Combine
  final_table <- rbind(female_table, male_table)
  
  # Save
  if (!dir.exists("results/tables")) dir.create("results/tables", recursive = TRUE)
  write.csv(final_table, "results/tables/table4.csv", row.names = FALSE)
  
  invisible(final_table)
}
