library(splines)
library(ggplot2)
library(dplyr)
library(ggpointdensity)
library(viridis)

# test_fig3 : MI vs IP standardized residual correlation plot
test_fig3 <- function() {
  
  # MI (increment) model with interaction
  mi_mod <- glm(
    INC ~ PL + INT * SEX,
    data = lobster,
    family = Gamma(link = "log")
  )
  
  # IP model with spline in PL
  ip_mod <- glm(
    INT ~ ns(PL, df = 3),
    data = lobster,
    family = Gamma(link = "log")
  )
  
  # Standardized residuals
  mi_res <- scale(residuals(mi_mod, type = "deviance"))[,1]
  ip_res <- scale(residuals(ip_mod, type = "deviance"))[,1]
  
  df <- data.frame(
    MI_res = mi_res,
    IP_res = ip_res,
    SEX    = lobster$SEX
  ) %>%
    filter(
      is.finite(MI_res),
      is.finite(IP_res),
      !is.na(SEX)
    )
  
  # Correlation by sex
 # Compute dynamic label positions within each facet
cor_df <- df %>%
  group_by(SEX) %>%
  summarise(
    r = round(cor(MI_res, IP_res), 2),
    x = quantile(MI_res, 0.80),   # 80% of x-range
    y = quantile(IP_res, 1.00)    # 90% of y-range
  )
 
  # Plot
  p <- ggplot(df, aes(MI_res, IP_res)) +
    geom_pointdensity(size = 2.0) +
    scale_color_viridis(option = "magma", direction = 1) +
    geom_smooth(method = "lm", se = FALSE,
                colour = "steelblue4", linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_text(
      data = cor_df,
      aes(x = x, y = y, label = paste0("r = ", r)),
      inherit.aes = FALSE,
      colour = "grey25",
      size = 4.5
    ) +
    facet_wrap(~ SEX, nrow = 1) +
    labs(
      title = "Standardized residuals for MI and IP",
      x = "Standardized MI residuals",
      y = "Standardized IP residuals"
    ) +
    theme_bw(base_size = 13) +
    theme(
      strip.background = element_rect(fill = "grey94", colour = "grey80"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  
  # Save
  if (!dir.exists("results/figures"))
    dir.create("results/figures", recursive = TRUE)
  
  ggsave("results/figures/fig3.png", p, width = 8, height = 5, dpi = 300)

  return(p)
}
