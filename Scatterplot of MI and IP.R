library(ggplot2)
library(dplyr)

# Step 1: Fit Gamma GLMs
mi_model <- glm(INC ~ PL + INT * SEX, family = Gamma(link = "log"), data = lobster)
ip_model <- glm(INT ~ PL + PL2, family = Gamma(link = "log"), data = lobster)

# Step 2: Extract residuals
mi_resid <- residuals(mi_model, type = "deviance")
ip_resid <- residuals(ip_model, type = "deviance")

# Step 3: Build dataframe (no link to `lobster`)
plot_data <- data.frame(
  MI_Residual = mi_resid,
  IP_Residual = ip_resid,
  SEX = lobster$SEX
)

# Step 4: Clean residuals + create factor variable safely
plot_data <- plot_data %>%
  filter(
    is.finite(MI_Residual),
    is.finite(IP_Residual),
    !is.na(SEX)
  ) %>%
  mutate(
    MI_Residual = (MI_Residual - mean(MI_Residual)) / sd(MI_Residual),
    IP_Residual = (IP_Residual - mean(IP_Residual)) / sd(IP_Residual),
    Sex = factor(ifelse(SEX == 1, "Female", "Male"), levels = c("Female", "Male"))
  )

# ✅ Step 5: Check that no non-finite values remain
cat("Total rows after cleaning:", nrow(plot_data), "\n")
cat("Any NA or NaN or Inf in MI_Residual?", any(!is.finite(plot_data$MI_Residual)), "\n")
cat("Any NA or NaN or Inf in IP_Residual?", any(!is.finite(plot_data$IP_Residual)), "\n")
cat("Any NA in Sex?", any(is.na(plot_data$Sex)), "\n")

# Step 6: Correlation text
cor_data <- plot_data %>%
  group_by(Sex) %>%
  summarize(r = cor(MI_Residual, IP_Residual), .groups = "drop") %>%
  mutate(
    x = 2.5,
    y = ifelse(Sex == "Female", 2, 1),
    label = paste0("r = ", round(r, 2))
  )

# ✅ Step 7: Plot using ONLY plot_data — no reference to `lobster`
ggplot(data = plot_data, mapping = aes(x = MI_Residual, y = IP_Residual, color = Sex)) +
  geom_point(alpha = 0.7, size = 2, na.rm = TRUE) +
  geom_smooth(method = "lm", se = FALSE, na.rm = TRUE) +
  geom_text(
    data = cor_data,
    mapping = aes(x = x, y = y, label = label, color = Sex),
    inherit.aes = FALSE,
    size = 5,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("blue", "red")) +
  labs(
    title = "Scatterplot of MI and IP residuals by sex",
    x = "Standardized MI Residuals",
    y = "Standardized IP Residuals",
    color = "Sex"
  ) +
  xlim(-3, 3) +
  ylim(-4, 3) +
  theme_minimal(base_size = 14)