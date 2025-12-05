library(ggplot2)

test_fig6 <- function(save_csv = TRUE) {
  
  ## 1. Validate dataset
  required_cols <- c("sex", "estimate", "lower", "upper",
                     "parameter_label", "Model")
  
  ## Optional: clean ordering for consistent facet layout
  ci_both$sex <- factor(ci_both$sex, levels = c("Female", "Male"))
  
  if (!is.factor(ci_both$parameter_label)) {
    ci_both$parameter_label <- factor(
      ci_both$parameter_label,
      levels = c("k", "L[infty]")
    )
  }
  
  ## 2. Build CI plot
  p_aoas <- ggplot(ci_both,
                   aes(x = sex, y = estimate,
                       shape = sex, linetype = sex)) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  width = 0.12, linewidth = 0.7, colour = "black") +
    geom_point(size = 2.8, colour = "black") +
    facet_grid(
      parameter_label ~ Model,
      scales   = "free_y",
      labeller = labeller(
        parameter_label = label_parsed,
        Model           = label_parsed
      )
    ) +
    labs(
      title = "Profile-likelihood 95% confidence intervals for k and L\u221e\nunder alternative random-effects distributions",
      x = "Sex",
      y = "Estimate"
    ) +
    
    # Female = open circle + solid CI; Male = filled circle + dashed CI
    scale_shape_manual(values = c("Female" = 1, "Male" = 16)) +
    scale_linetype_manual(values = c("Female" = "solid", "Male" = "dashed")) +
    
    theme_bw(base_size = 12) +
    theme(
      plot.title   = element_text(size = 12, hjust = 0.5),
      strip.text   = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "grey95"),
      axis.title   = element_text(size = 11),
      axis.text    = element_text(size = 10),
      legend.title = element_blank(),
      legend.position = "top"
    )
  
  ## 3. Save plot
  if (!dir.exists("results/figures")) {
    dir.create("results/figures", recursive = TRUE)
  }
  
  ggsave("results/figures/fig6.png",
         plot = p_aoas,
         width = 6, height = 4, dpi = 300)
  
  ## 4. (Optional) Save CI dataset for reproducibility
  if (save_csv) {
    if (!dir.exists("results/tables")) {
      dir.create("results/tables", recursive = TRUE)
    }
    write.csv(ci_both,
              "results/tables/fig6_ci_dataset.csv",
              row.names = FALSE)
  }
  
  ## Return plot object invisibly
  invisible(p_aoas)
}
