# Test functions orchestrator
# This script runs basic checks to ensure reproducibility

# ✅ Check that data file exists
stopifnot(file.exists("data/lobster.csv"))

# ✅ Check that core analysis scripts can be sourced
source("scripts/fig3.R")
source("scripts/fig4.R")
source("scripts/fig6.R")
source("scripts/fig7.R")

# ✅ Check that table scripts run without error
source("scripts/tables/table2_bf.R")
source("scripts/tables/table2_igf.R")
source("scripts/tables/table2_gf.R")
source("scripts/tables/table2_bg_model.R")
source("scripts/tables/table2_bl_model.R")
source("scripts/tables/table2_gg_model.R")
source("scripts/tables/table2_gl_model.R")
source("scripts/tables/table3_bg_model.R")
source("scripts/tables/table3_bl_model.R")
source("scripts/tables/table4.R")
source("scripts/tables/table5.R")
source("scripts/tables/table6_independent_model.R")
source("scripts/tables/table6_joint_models.R")

# ✅ Check that plot scripts run without error
source("scripts/plots/fig1.R")
source("scripts/plots/fig2.R")
source("scripts/plots/fig5.R")

cat("All tests passed successfully!\n")
