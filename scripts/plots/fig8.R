library(ggplot2)
library(gridExtra)

## test_fig9: smoke test for MI residual diagnostics
test_fig8 <- function() {
  
  # Read data
  df <- lobster
  
  # Compute proportion of growth relative to next size
  epsilon <- 1e-5
  df$ratio <- with(df, INC / (CL + INC))
  df$ratio[df$ratio < epsilon] <- epsilon
  df$ratio[df$ratio > 1 - epsilon] <- 1 - epsilon
  
  # Logistic transformation and linear regression on premoult length
  logit_ratio <- log(df$ratio / (1 - df$ratio))
  model <- lm(logit_ratio ~ CL, data = df)
  
  # Predicted ratio and standardized residuals
  df$pred_ratio  <- plogis(predict(model))
  df$std_resid   <- (df$ratio - df$pred_ratio) /
    sqrt(df$pred_ratio * (1 - df$pred_ratio))
  df$abs_resid   <- abs(df$std_resid)
  df$scale_resid <- sqrt(df$abs_resid)
  df$index       <- seq_len(nrow(df))
  
  # Define point aesthetics
  point_opts <- list(color = "blue", alpha = 0.7, size = 1.2)
  
  # 1. Residuals vs fitted
  p1 <- ggplot(df, aes(x = pred_ratio, y = std_resid)) +
    do.call(geom_point, point_opts) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.8) +
    labs(title = "MI: Residuals vs fitted",
         x = "Fitted ratio",
         y = "Standardized residuals") +
    theme_minimal()
  
  # 2. Scale–location plot
  p2 <- ggplot(df, aes(x = pred_ratio, y = scale_resid)) +
    do.call(geom_point, point_opts) +
    geom_smooth(method = "loess", se = FALSE, colour = "red", size = 1.5) +
    labs(title = "MI: Scale–location",
         x = "Fitted ratio",
         y = "sqrt(|residuals|)") +
    theme_minimal()
  
  # 3. Normal Q–Q plot
  p3 <- ggplot(df, aes(sample = std_resid)) +
    stat_qq(color = point_opts$color, alpha = point_opts$alpha, size = point_opts$size) +
    stat_qq_line(colour = "red", size = 1.5) +
    labs(title = "MI: Normal Q–Q plot") +
    theme_minimal()
  
  # 4. Influence proxy plot
  p4 <- ggplot(df, aes(x = index, y = abs_resid)) +
    do.call(geom_point, point_opts) +
    geom_hline(yintercept = 3, linetype = "dashed", size = 0.8) +
    labs(title = "MI: Influence proxy",
         x = "Observation index",
         y = "|Standardized residual|") +
    theme_minimal()
  
  # Save panel
  if (!dir.exists("results/figures")) dir.create("results/figures", recursive = TRUE)
  png("results/figures/fig8.png", width = 1200, height = 900, res = 150)
  grid.arrange(p1, p2, p3, p4, ncol = 2)
  dev.off()
  
  invisible(list(model = model, plots = list(p1, p2, p3, p4)))
}
