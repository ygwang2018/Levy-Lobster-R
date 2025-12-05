library(ggplot2)
library(dplyr)

## ---------------------------------------------------------
## test_fig4: smoke test for lobster growth trajectories
## ---------------------------------------------------------
test_fig4 <- function() {
  cat("Running test for fig4...\n")
  
  # 1. Convert date
  lobster$Date <- as.Date(lobster$Date, "%d/%m/%Y")
  
  # 2. Factors
  lobster$SEX     <- factor(lobster$SEX, levels = c("1", "2"), labels = c("Female", "Male"))
  lobster$LOBSTER <- factor(lobster$LOBSTER)
  
  # 3. Assign colors per lobster (consistent within facet)
  n_lobsters <- length(levels(lobster$LOBSTER))
  lobster_colors <- scales::hue_pal()(n_lobsters)
  names(lobster_colors) <- levels(lobster$LOBSTER)
  
  # 4. Plot
  p <- ggplot(lobster, aes(x = Date, y = CL, group = LOBSTER, color = LOBSTER)) +
    geom_step(linewidth = 0.7, alpha = 0.90) +
    facet_wrap(~SEX, nrow = 1, scales = "free_x",
               strip.position = "top") +
    scale_color_manual(values = lobster_colors) +
    labs(
      x = "Time (year)",
      y = "Carapace Length (mm)"
    ) +
    scale_x_date(date_labels = "%Y") +
    theme_bw(base_size = 14) +
    theme(
      strip.background = element_rect(fill = "peachpuff", color = NA),
      strip.text = element_text(face = "bold", size = 15),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey85")
    )
  
  # 5. Save to results/figures
  if (!dir.exists("results/figures")) dir.create("results/figures", recursive = TRUE)
  ggsave("results/figures/fig4.png", plot = p, width = 8, height = 5, dpi = 300)
  
  cat("fig4.png saved to results/figures.\n")
  
  invisible(p)
}
