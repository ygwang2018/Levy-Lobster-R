# MI ~ PL + INT * SEX
mi_model <- glm(
    INC ~ PL + INT * SEX,
    data = lobster,
    family = Gamma(link = "log")
)

# IP ~ PL + PL^2
ip_model <- glm(
    INT ~ PL + I(PL^2),
    data = lobster,
    family = Gamma(link = "log")
)
library(splines)

ip_model <- glm(
    INT ~ ns(PL, df = 3),
    data = lobster,
    family = Gamma(link = "log")
)
library(ggplot2)
library(dplyr)
library(ggpointdensity)
library(viridis)

# 1. Compute standardized residuals for MI and IP
mi_res  <- residuals(mi_model, type = "deviance")
ip_res  <- residuals(ip_model, type = "deviance")

plot_data <- data.frame(
    MI_res = scale(mi_res)[,1],
    IP_res = scale(ip_res)[,1],
    SEX = lobster$SEX
) %>%
    filter(
        is.finite(MI_res),
        is.finite(IP_res),
        !is.na(SEX)
    )

# 2. Compute Pearson correlations per sex
cor_data <- plot_data %>%
    group_by(SEX) %>%
    summarise(r = round(cor(MI_res, IP_res), 2)) %>%
    mutate(
        x = 2.5,   # label x-position
        y = 3.7    # label y-position
    )

# Assign plot to an object
p <- ggplot(plot_data, aes(x = MI_res, y = IP_res)) +
  geom_pointdensity(size = 2.2) +
  scale_color_viridis(
    option = "magma",
    direction = 1,
    name = "Point density"
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    colour = "steelblue4",
    linewidth = 0.9
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey65") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey65") +
  geom_text(
    data = cor_data,
    aes(x = x, y = y, label = paste0("r = ", r)),
    inherit.aes = FALSE,
    colour = "grey25",
    size = 4.5,
    hjust = 1
  ) +
  facet_wrap(~SEX, nrow = 1, strip.position = "top") +
  labs(
    title = "Standardized residuals for MI and IP",
    x = "Standardized MI residuals",
    y = "Standardized IP residuals"
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.background = element_rect(fill = "grey95", colour = "grey85"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    panel.grid.major = element_line(colour = "grey90"),
    panel.grid.minor = element_blank()
  )

# Save to results/figures
ggsave("results/figures/mi_ip_residuals.png", plot = p,
       width = 8, height = 5, dpi = 300)
