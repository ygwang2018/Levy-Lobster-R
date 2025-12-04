# Main analysis orchestrator
# This script sources all table and analysis scripts

source("scripts/lobster_growth.R")
source("scripts/mcmc_simulation.R")
source("scripts/profile_likelihood_ci.R")
source("scripts/obs_vs_predicted_ip.R")

# Tables
source("scripts/tables/table2_bg_model.R")
source("scripts/tables/table2_gl_model.R")
source("scripts/tables/table2_gamma_linf_fixed.R")
source("scripts/tables/table2_beta_linf_fixed.R")
source("scripts/tables/table6_joint_models.R")
# Main analysis script placeholder
