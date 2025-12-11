test_fig1 <- function() {

  ## Prepare data
  df <- lobster
  df$Sex <- factor(df$SEX, levels = c(1, 2), labels = c("Female", "Male"))

  ## Colours used in manuscript
  cols <- c(Female = "#1f78b4",
            Male   = "#e31a1c")

  ## Plot intermoult period distribution by sex
  p <- ggplot(df, aes(x = INT, fill = Sex)) +
    geom_histogram(
      binwidth = 5,
      boundary = 0,
      colour   = "black",
      alpha    = 0.85
    ) +
    facet_grid(Sex ~ ., scales = "fixed") +
    scale_fill_manual(values = cols) +
    labs(
      title = "Intermoult period distribution by sex",
      x     = "Intermoult period (days)",
      y     = "Count"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title        = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor  = element_blank(),
      panel.grid.major.x= element_blank(),
      axis.line         = element_line(colour = "black"),
      strip.text        = element_text(size = 13, face = "bold"),
      legend.position   = "none",
      panel.spacing     = unit(1.2, "lines")
    )

  ## Save figure
  if (!dir.exists("results/figures"))
    dir.create("results/figures", recursive = TRUE)

  ggsave("results/figures/fig1.png", p,
         width = 6, height = 4, dpi = 300)

  return(p)
}
