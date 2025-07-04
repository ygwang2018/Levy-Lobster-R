# Plot
ggplot(ci_data, aes(x=sex, y=estimate, color=sex)) +
geom_point(size=4) +
geom_errorbar(aes(ymin=lower, ymax=upper), width=0.15, linewidth=1) +
facet_wrap(~parameter_label, scales="free_y", labeller = label_parsed) +
labs(
title = "Profile Likelihood Confidence Intervals by Sex",
x = "Sex",
y = "Estimate"
) +
scale_color_manual(values=c("#D55E00", "#0072B2")) +
theme_minimal(base_size = 14) +
theme(
strip.text = element_text(size=16, face="bold"),
legend.title = element_blank(),
legend.position = "top",
axis.text = element_text(size=12)
)
