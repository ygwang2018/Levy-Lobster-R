## 1. Correct SEX factor
lobster$SEX <- factor(lobster$SEX,
                      levels = c("1", "2"),
                      labels = c("Female", "Male"))

## 2. Fit model
glm_int <- glm(
  INT ~ ns(PL, df = 3) + ns(NINT, df = 3),
  data   = lobster,
  family = Gamma(link = "log")
)

## 3. Predictions
lobster$Predicted_INT <- predict(glm_int, type = "response")

## 4. Remove values outside plotting scale
lobster_clean <- lobster %>%
  filter(
    INT >= 0, INT <= 500,
    Predicted_INT >= 0, Predicted_INT <= 500
  )

## 5. Plot
ggplot(lobster_clean, aes(x = INT, y = Predicted_INT, colour = SEX)) +
  geom_point(alpha = 0.7, size = 2.3) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed",
              linewidth = 1,
              colour = "grey40") +
  facet_wrap(~SEX, nrow = 1, strip.position = "top") +
  scale_colour_manual(values = c("Female" = "#F28E2B",
                                 "Male"   = "#4E79A7")) +
  scale_x_continuous(limits = c(0, 500)) +
  scale_y_continuous(limits = c(0, 500)) +
  theme_bw()
