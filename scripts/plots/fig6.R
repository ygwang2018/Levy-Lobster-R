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
  scale_shape_manual(values = c("Female" = 1,   # open circle
                                "Male"   = 16)) +  # filled circle
  scale_linetype_manual(values = c("Female" = "solid",
                                   "Male"   = "dashed")) +
  theme_bw(base_size = 12) +
  theme(
    plot.title   = element_text(size = 12, hjust = 0.5),
    strip.text   = element_text(size = 11, face = "bold"),
    axis.title   = element_text(size = 11),
    axis.text    = element_text(size = 10),
    legend.title = element_blank(),
    legend.position = "top"
  )

print(p_aoas)

# Save plot to results/figures
ggsave("results/figures/fig6.png", plot = p_aoas, width = 6, height = 4, dpi = 300)

