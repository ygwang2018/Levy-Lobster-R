library(ggplot2)

# Load raw data
data <- lobster

# Convert SEX
data$Sex <- factor(data$SEX, levels = c(1, 2), labels = c("Female", "Male"))

# Plot
ggplot(data, aes(x = INT, color = Sex, linetype = Sex)) +
  geom_density(size = 1) +
  labs(
    title = "Intermoult period distribution",
    x = "Time (days)",
    y = "Density"
  ) +
  scale_color_manual(values = c("Female" = "black", "Male" = "red")) +
  scale_linetype_manual(values = c("Female" = "solid", "Male" = "dashed")) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "top"
  )
