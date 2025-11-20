ggplot(plot_data, aes(x = MI_Residual, y = IP_Residual)) +
  
  geom_pointdensity(size = 2.2) +
  
  scale_color_viridis(
    option = "magma",
    direction = 1,    # ✔ purple → pink → orange → yellow
    name = "Point density"
  ) +
  
  geom_smooth(
    method = "lm",
    se = FALSE,
    colour = "steelblue4",
    size = 0.9
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  
  geom_text(
    data = cor_data,
    aes(x = x, y = y, label = label),
    hjust = 1.2, vjust = 1.4,
    colour = "grey20",
    size = 4.5,
    inherit.aes = FALSE
  ) +
  
  facet_wrap(~ Sex, nrow = 1) +
  
  labs(
    title = "Standardized residuals for MI and IP",
    x = "Standardized MI residuals",
    y = "Standardized IP residuals"
  ) +
  
  theme_bw(base_size = 15) +
  theme(
    panel.grid.major = element_line(colour = "grey90"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = "grey80"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )
