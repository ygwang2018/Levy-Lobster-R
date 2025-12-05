library(ggplot2)

test_fig6 <- function() {
  cat("Running test for fig6...\n")
  
  # Check dataset exists
  if (!exists("ci_both")) stop("Dataset 'ci_both' not found.")
  if (!all(c("sex","estimate","lower","upper","parameter_label","Model") %in% names(ci_both))) {
    stop("Missing required columns in 'ci_both'.")
  }
  
  # Build plot
  p_aoas <- ggplot(ci_both,
                   aes(x = sex, y = estimate,
                       shape = sex, linetype = sex)) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  width = 0.12, linewidth = 0.7, colour = "black") +
    geom_point(size = 2.8, colour = "black") +
    facet_grid(parameter_label ~ Model,
               scales  = "free_y",
               labeller = labeller(
                 parameter_label = label_parsed,
                 Model           = label_parsed
               )) +
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
      axis.title   = element_text(size = 11),
      axis.text    = element_text(size = 10),
      legend.title = element_blank(),
      legend.position = "top"
    )
  
  # Save plot
  if (!dir.exists("results/figures")) dir.create("results/figures", recursive = TRUE)
  ggsave("results/figures/fig6.png", plot = p_aoas, width = 6, height = 4, dpi = 300)
  
  invisible(p_aoas)
}
