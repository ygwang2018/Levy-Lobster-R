library(ggplot2)
library(dplyr)


ci_both <- bind_rows(

  # -------------------- LOGNORMAL MODEL --------------------
  data.frame(
    Model = "Lognormal~L[infty]",
    parameter_label = "k",
    sex       = c("Female", "Male"),
    estimate  = c(0.277, 0.304),
    lower     = c(0.218, 0.254),
    upper     = c(0.335, 0.353)
  ),

  data.frame(
    Model = "Lognormal~L[infty]",
    parameter_label = "L[infty]",
    sex       = c("Female", "Male"),
    estimate  = c(185.32, 201.48),
    lower     = c(163.75, 182.02),
    upper     = c(206.88, 220.94)
  ),

  # -------------------- GAMMA MODEL ------------------------
  data.frame(
    Model = "Gamma~L[infty]",
    parameter_label = "k",
    sex       = c("Female", "Male"),
    estimate  = c(0.267, 0.337),
    lower     = c(0.246, 0.313),
    upper     = c(0.288, 0.360)
  ),

  data.frame(
    Model = "Gamma~L[infty]",
    parameter_label = "L[infty]",
    sex       = c("Female", "Male"),
    estimate  = c(189.80, 190.11),
    lower     = c(186.39, 185.84),
    upper     = c(193.19, 195.28)
  )
)

# Ordering of factors
ci_both$sex <- factor(ci_both$sex, levels = c("Female", "Male"))
ci_both$parameter_label <- factor(ci_both$parameter_label, levels = c("k", "L[infty]"))


test_fig6 <- function(save_csv = TRUE) {

  # Build plot
  p_aoas <- ggplot(ci_both,
                   aes(x = sex, y = estimate,
                       shape = sex, linetype = sex)) +

    geom_errorbar(aes(ymin = lower, ymax = upper),
                  width = 0.12, linewidth = 0.7, colour = "black") +

    geom_point(size = 3, colour = "black") +

    facet_grid(parameter_label ~ Model,
               scales   = "free_y",
               labeller = labeller(
                 parameter_label = label_parsed,
                 Model           = label_parsed
               )) +

    scale_shape_manual(values = c("Female" = 1, "Male" = 16)) +
    scale_linetype_manual(values = c("Female" = "solid", "Male" = "dashed")) +

    labs(
      x = "Sex",
      y = "Estimate",
      title = NULL
    ) +

    theme_bw(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_rect(fill = "grey92"),
      legend.position = "top",
      legend.title = element_blank()
    )

  # Save figure
  dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
  ggsave("results/figures/fig5.png", plot = p_aoas,
         width = 7, height = 5, dpi = 300)

  # Save CI dataset
  if (save_csv) {
    dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
    write.csv(ci_both, "results/tables/fig5_ci_dataset.csv", row.names = FALSE)
  }

  p_aoas
}


test_fig6()
