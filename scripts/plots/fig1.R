library(ggplot2)

df <- lobster

# Convert SEX to factor safely
df$Sex <- factor(df$SEX)
if (length(levels(df$Sex)) == 2) {
  levels(df$Sex) <- c("Female", "Male")
}

# Choose elegant, soft colours
cols <- c("Female" = "#1f78b4",   # soft blue
          "Male"   = "#e31a1c")   # deep muted red


test_fig1 <- function() {
  cat("Running test_function for fig2...\n")
  
  if (!exists("lobster")) stop("Dataset 'lobster' not found.")
  if (!all(c("SEX","INT") %in% names(lobster))) stop("Missing required columns.")
  
  # Re-run the plotting code quickly
  df <- lobster
  df$Sex <- factor(df$SEX)
  if (length(levels(df$Sex)) == 2) levels(df$Sex) <- c("Female","Male")
  
  cols <- c("Female"="#1f78b4","Male"="#e31a1c")
  
  p <- ggplot(df, aes(x=INT, fill=Sex)) +
    geom_histogram(binwidth=5, boundary=0, colour="black", alpha=0.85, na.rm=TRUE) +
    facet_grid(Sex ~ ., scales="fixed") +
    scale_fill_manual(values=cols) +
    labs(title="Intermoult period distribution by sex",
         x="Intermoult period (days)", y="Count") +
    theme_minimal(base_size=14)
  
  if (!dir.exists("results/figures")) dir.create("results/figures", recursive=TRUE)
  ggsave("results/figures/fig2.png", plot=p, width=6, height=4, dpi=300)
}



# Plot: vertical stacked histograms with colour + elegance
p<- ggplot(df, aes(x = INT, fill = Sex)) +
  geom_histogram(
    binwidth = 5,
    boundary = 0,
    colour   = "black",
    alpha    = 0.85,
    na.rm    = TRUE
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
    plot.title = element_text(
      hjust = 0.5,
      face  = "bold"
    ),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line          = element_line(colour = "black"),
    strip.text         = element_text(size = 13, face = "bold"),
    legend.position    = "none",     # clean: sex appears as facet labels
    panel.spacing      = unit(1.2, "lines")
  )
# Save plot to results/figures
ggsave("results/figures/fig1.png", plot = p, width = 6, height = 4, dpi = 300)
