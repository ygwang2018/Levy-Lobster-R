library(dplyr)
library(splines)

# Standardize factor levels
lobster$SEX <- factor(lobster$SEX)

# Optional: relabel levels if needed
levels(lobster$SEX) <- c("Female", "Male")  # Match your intended labels

# Fit model
glm_int <-glm(INT ~ ns(PL, df = 3) + ns(NINT, df = 3),
    data = lobster, family = Gamma(link = "log"))

# Predict
lobster$Predicted_INT <- predict(glm_int, type = "response")

# Plot
ggplot(lobster, aes(x = INT, y = Predicted_INT, color = SEX)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~SEX) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Female" = "red", "Male" = "cyan")) +
  coord_cartesian(ylim = c(0, 500)) +
  labs(
    title = "Observed vs Predicted INT",
    x = "Observed INT",
    y = "Predicted INT",
    color = "SEX"
  ) +
  theme_minimal()
