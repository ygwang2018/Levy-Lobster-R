test_fig4 <- function() {

 
  # convert date
  lobster$Date <- as.Date(lobster$Date, "%d/%m/%Y")

  # factors
  lobster$SEX <- factor(lobster$SEX, levels = c(1, 2),
                        labels = c("Female", "Male"))
  lobster$LOBSTER <- factor(lobster$LOBSTER)

  # colours for individuals
  n_id <- length(levels(lobster$LOBSTER))
  id_cols <- scales::hue_pal()(n_id)
  names(id_cols) <- levels(lobster$LOBSTER)

  # plot
  p <- ggplot(lobster,
              aes(x = Date, y = CL,
                  group = LOBSTER,
                  colour = LOBSTER)) +
    geom_step(linewidth = 0.7, alpha = 0.9) +
    facet_wrap(~ SEX, nrow = 1, scales = "free_x") +
    scale_colour_manual(values = id_cols) +
    scale_x_date(date_labels = "%Y") +
    labs(x = "Time (year)", y = "Carapace Length (mm)") +
    theme_bw(base_size = 14) +
    theme(strip.background = element_rect(fill = "peachpuff", colour = NA),
          strip.text = element_text(face = "bold", size = 15),
          legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey85"))

  # save
  if (!dir.exists("results/figures"))
    dir.create("results/figures", recursive = TRUE)

  ggsave("results/figures/fig4.png",
         plot = p, width = 8, height = 5, dpi = 300)

  invisible(p)
  print(p)
}


test_fig4()
