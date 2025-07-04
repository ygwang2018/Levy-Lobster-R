library(ggplot2)

# Load your data
data <- lobster

# Convert SEX
data$Sex <- factor(data$SEX, levels = c(1, 2), labels = c("Female", "Male"))

# Remove missing or invalid values
data_clean <- subset(data, !is.na(INC))

# Plot
ggplot(data_clean, aes(x = INC, color = Sex, linetype = Sex)) +
  geom_density(size = 1) +
  labs(
    title = "Increment (INC) distribution",
    x = "Increment (INC)",
    y = "Density"
  ) +
  scale_color_manual(values = c("Female" = "black", "Male" = "red")) +
  scale_linetype_manual(values = c("Female" = "solid", "Male" = "dashed")) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "top"
  )
