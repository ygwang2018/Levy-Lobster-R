# Test functions orchestrator
# This script runs basic checks to ensure reproducibility

# ✅ Check that data file exists
stopifnot(file.exists("data/lobster.csv"))

# ✅ Check that core analysis scripts can be sourced
source("scripts/fig1.R")
source("scripts/fig3.R")
source("scripts/fig4.R")
source("scripts/fig5.R")

# ✅ Check that table scripts run without error
source("scripts/tables/table2_bf.R")
source("scripts/tables/table2_igf.R")
source("scripts/tables/table2_gf.R")
source("scripts/tables/table2_bg_model.R")
source("scripts/tables/table2_bl_model.R")
source("scripts/tables/table2_gg_model.R")
source("scripts/tables/table2_gl_model.R")
source("scripts/tables/supplementary_bg_model.R")
source("scripts/tables/supplementary_bl_model.R")
source("scripts/tables/table3.R")
source("scripts/tables/table4.R")
source("scripts/tables/table5_independent_model.R")
source("scripts/tables/table5_joint_models.R")

cat("All tests passed successfully!\n")
