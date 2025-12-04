# Test functions orchestrator
# This script runs basic checks to ensure reproducibility

# ✅ Check that data file exists
stopifnot(file.exists("data/lobster.csv"))

# ✅ Check that core analysis scripts can be sourced
source("scripts/fig4.R")
source("scripts/fig5.R")
source("scripts/fig8.R")
source("scripts/fig10.R")

# ✅ Check that table scripts run without error
source("scripts/tables/table2_bg_model.R")
source("scripts/tables/table2_gl_model.R")
source("scripts/tables/table2_gamma_linf_fixed.R")
source("scripts/tables/table2_beta_linf_fixed.R")
source("scripts/tables/table6_joint_models.R")

# ✅ Check that plot scripts run without error
source("scripts/plots/fig2.R")
source("scripts/plots/fig3.R")
source("scripts/plots/fig4.R")
source("scripts/plots/fig5.R")
source("scripts/plots/fig6.R")
source("scripts/plots/fig7.R")
source("scripts/plots/fig8.R")
source("scripts/plots/fig9.R")
source("scripts/plots/fig10.R")

cat("All tests passed successfully!\n")
