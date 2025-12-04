library(ggplot2)
library(dplyr)

# Custom elegant colours (journal-safe)
sex_cols <- c("Female" = "#3AAFA9",   # soft teal
              "Male"   = "#8E44AD")   # elegant plum

## ---------------------------------------------------------
## test_function: smoke test for fig3 plotting script
## ---------------------------------------------------------
test_function <- function() {
  cat("Running test_function for fig3...\n")
  
  # Check dataset exists
  if (!exists("renew1")) stop("Dataset 'renew1' not found.")
  if (!all(c("SEX","INC") %in% names(renew1))) stop("Missing required columns.")
  
  # Prepare data
  df <- renew1 %>%
    mutate(SEX = factor(SEX, levels = c(1, 2), labels = c("Female", "Male")))
  
  # Build plot
  p <- ggplot(df, aes(x = INC, fill = SEX)) +
    geom_histogram(bins = 30, colour = "white", linewidth = 0.25, alpha = 0.9) +
    facet_wrap(~ SEX, nrow = 2, strip.position = "top") +
    scale_fill_manual(values = sex_cols) +
    labs(title = "Moult increment distribution",
         x = "Increment (mm)", y = "Count") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      strip.background = element_rect(fill = "grey94", colour = "grey80"),
      strip.text = element_text(face = "bold", size = 13),
      panel.grid.major = element_line(colour = "grey88", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(colour = "black"),
      panel.spacing = unit(1.3, "lines"),
      legend.position = "none"
    )
  
  # Save plot
  if (!dir.exists("results/figures")) dir.create("results/figures", recursive = TRUE)
  ggsave("results/figures/fig3.png", plot = p, width = 6, height = 4, dpi = 300)
  
  invisible(p)
}
