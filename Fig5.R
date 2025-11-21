library(lattice)

# 1. Convert date
lobster$Date <- as.Date(lobster$Date, "%d/%m/%Y")

# 2. Assign factor labels
lobster$SEX <- factor(lobster$SEX, levels = c("1", "2"), labels = c("Female", "Male"))
lobster$LOBSTER <- factor(lobster$LOBSTER)

# 3. Assign colors to each lobster
n_lobsters <- length(levels(lobster$LOBSTER))
lobster_colors <- rainbow(n_lobsters)
names(lobster_colors) <- levels(lobster$LOBSTER)

# 4. Plot
xyplot(
  CL ~ Date | SEX,
  data = lobster,
  type = "s",
  groups = LOBSTER,
  layout = c(2, 1),
  ylim = c(0, 170),
  xlab = "Time (year)",
  ylab = "Carapace Length (mm)",
  strip = strip.custom(bg = "peachpuff", factor.levels = c("Female", "Male")),
  scales = list(
    x = list(format = "%Y", relation = "same"),
    y = list(relation = "same")
  ),
  panel = function(x, y, groups, subscripts, ...) {
    gp <- groups[subscripts]
    for (lob_id in unique(gp)) {
      idx <- which(gp == lob_id)
      col <- lobster_colors[as.character(lob_id)]
      panel.lines(x[idx], y[idx], type = "s", col = col)
    }
  }
)
