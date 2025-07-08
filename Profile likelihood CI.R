# ----------------------------------------
# Load required libraries
# ----------------------------------------
library(ggplot2)
library(dplyr)
library(numDeriv)
dat<-lobster
# ----------------------------------------
# Define functions
# ----------------------------------------
ff <- function(Linf, k, zeta, miu, sigma, Xi, Ti, Li) {
  gLinf <- 0
  ni <- length(Xi)
  for (j in 1:ni) {
    a <- (1 - exp(-k * Ti[j])) * (zeta - 1)
    b <- exp(-k * Ti[j]) * (zeta - 1)
    x <- Xi[j] / (Linf - Li[j])
    bb <- lgamma(a + b) + (a - 1) * log(x) + (b - 1) * log(1 - x) - lgamma(a) - lgamma(b) -
      log(sigma) - 0.5 * log(2 * pi) - ((log(Linf) - miu)^2) / (2 * sigma^2) - log(Linf - Li[j])
    gLinf <- gLinf + bb
  }
  return(gLinf)
}

gg <- function(Linf, k, zeta, miu, sigma, Xi, Ti, Li, MinLinf, MaxLinf) {
  res <- ff(Linf, k, zeta, miu, sigma, Xi, Ti, Li)
  Mi <- optimize(ff, c(MinLinf, MaxLinf), k = k, zeta = zeta, miu = miu, sigma = sigma,
                 Xi = Xi, Ti = Ti, Li = Li, maximum = TRUE)$objective
  return(exp(res - Mi))
}

LL <- function(theta, dat, MinLinf, MaxLinf) {
  k <- theta[1]
  zeta <- theta[2]
  miu <- theta[3]
  sigma <- theta[4]
  all.LL <- 0
  for (iid in unique(dat$LOBSTER)) {
    dati <- dat[dat$LOBSTER == iid, ]
    Xi <- dati$INC
    Li <- dati$PL
    Ti <- dati$INT / 365.25
    Mi <- optimize(ff, c(MinLinf, MaxLinf), k = k, zeta = zeta, miu = miu, sigma = sigma,
                   Xi = Xi, Ti = Ti, Li = Li, maximum = TRUE)$objective
    res.int <- integrate(gg, lower = MinLinf, upper = MaxLinf, k = k, zeta = zeta, miu = miu,
                         sigma = sigma, Xi = Xi, Ti = Ti, Li = Li,
                         MinLinf = MinLinf, MaxLinf = MaxLinf)$value
    res.int <- log(res.int)
    all.LL <- all.LL + res.int + Mi
  }
  return(-all.LL)
}

# ----------------------------------------
# Common settings
# ----------------------------------------
MinLinf <- 160
MaxLinf <- 200
init_theta <- c(0.1, 4, 1, 2.88)
lb <- c(0.01, 0.1, 0.1, 0.1)
ub <- c(1, 100, 10, 10)

# ----------------------------------------
# Female
# ----------------------------------------
dat_f <- dat[dat$SEX == 1, ]

res_f <- optim(init_theta, LL, dat = dat_f,
               MinLinf = MinLinf, MaxLinf = MaxLinf,
               lower = lb, upper = ub,
               method = "L-BFGS-B",
               hessian = FALSE,
               control = list(maxit = 1000, trace = 1))

# Numerical Hessian
hess_f <- hessian(func = LL, x = res_f$par, dat = dat_f, MinLinf = MinLinf, MaxLinf = MaxLinf)
cov_mat_f <- tryCatch({ solve(hess_f) }, error = function(e) diag(length(res_f$par)) * NA)
se_f <- sqrt(diag(cov_mat_f))

# k
estimate_k_f <- res_f$par[1]
se_k_f <- se_f[1]
lower_k_f <- estimate_k_f - 1.96 * se_k_f
upper_k_f <- estimate_k_f + 1.96 * se_k_f

# Linf
mu_f <- res_f$par[3]
sigma_f <- res_f$par[4]
Linf_f <- exp(mu_f + (sigma_f^2)/2)

var_mu_f <- cov_mat_f[3, 3]
var_sigma_f <- cov_mat_f[4, 4]
cov_mu_sigma_f <- cov_mat_f[3, 4]

dLinf_dmu_f <- Linf_f
dLinf_dsigma_f <- Linf_f * sigma_f

var_Linf_f <- dLinf_dmu_f^2 * var_mu_f + dLinf_dsigma_f^2 * var_sigma_f + 2 * dLinf_dmu_f * dLinf_dsigma_f * cov_mu_sigma_f
se_Linf_f <- sqrt(var_Linf_f)

lower_Linf_f <- Linf_f - 1.96 * se_Linf_f
upper_Linf_f <- Linf_f + 1.96 * se_Linf_f

# ----------------------------------------
# Male
# ----------------------------------------
dat_m <- dat[dat$SEX == 2, ]

res_m <- optim(init_theta, LL, dat = dat_m,
               MinLinf = MinLinf, MaxLinf = MaxLinf,
               lower = lb, upper = ub,
               method = "L-BFGS-B",
               hessian = FALSE,
               control = list(maxit = 1000, trace = 1))

# Numerical Hessian
hess_m <- hessian(func = LL, x = res_m$par, dat = dat_m, MinLinf = MinLinf, MaxLinf = MaxLinf)
cov_mat_m <- tryCatch({ solve(hess_m) }, error = function(e) diag(length(res_m$par)) * NA)
se_m <- sqrt(diag(cov_mat_m))

# k
estimate_k_m <- res_m$par[1]
se_k_m <- se_m[1]
lower_k_m <- estimate_k_m - 1.96 * se_k_m
upper_k_m <- estimate_k_m + 1.96 * se_k_m

# Linf
mu_m <- res_m$par[3]
sigma_m <- res_m$par[4]
Linf_m <- exp(mu_m + (sigma_m^2)/2)

var_mu_m <- cov_mat_m[3, 3]
var_sigma_m <- cov_mat_m[4, 4]
cov_mu_sigma_m <- cov_mat_m[3, 4]

dLinf_dmu_m <- Linf_m
dLinf_dsigma_m <- Linf_m * sigma_m

var_Linf_m <- dLinf_dmu_m^2 * var_mu_m + dLinf_dsigma_m^2 * var_sigma_m + 2 * dLinf_dmu_m * dLinf_dsigma_m * cov_mu_sigma_m
se_Linf_m <- sqrt(var_Linf_m)

lower_Linf_m <- Linf_m - 1.96 * se_Linf_m
upper_Linf_m <- Linf_m + 1.96 * se_Linf_m

# ----------------------------------------
# Create summary data frame
# ----------------------------------------
ci_data <- data.frame(
  sex = c("Female", "Female", "Male", "Male"),
  parameter = c("k", "Linf", "k", "Linf"),
  estimate = c(estimate_k_f, Linf_f, estimate_k_m, Linf_m),
  lower = c(lower_k_f, lower_Linf_f, lower_k_m, lower_Linf_m),
  upper = c(upper_k_f, upper_Linf_f, upper_k_m, upper_Linf_m)
)

# Add expression labels
ci_data <- ci_data %>%
  mutate(parameter_label = case_when(
    parameter == "k" ~ "italic(k)",
    parameter == "Linf" ~ "L[infinity]",
    TRUE ~ parameter
  ))

# ----------------------------------------
# Plot
# ----------------------------------------
ggplot(ci_data, aes(x = sex, y = estimate, color = sex)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15, linewidth = 1) +
  facet_wrap(~parameter_label, scales = "free_y", labeller = label_parsed) +
  labs(
    title = "Profile Likelihood 95% Confidence Intervals by Sex",
    x = "Sex",
    y = "Estimate"
  ) +
  scale_color_manual(values = c("#D55E00", "#0072B2")) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    legend.title = element_blank(),
    legend.position = "top",
    axis.text = element_text(size = 12)
  )

# ----------------------------------------
# Display numeric results for female
# ----------------------------------------
cat("\n========================\n")
cat("Female estimates:\n")
cat("k estimate    =", round(estimate_k_f, 4), "\n")
cat("k 95% CI      =", round(lower_k_f, 4), "-", round(upper_k_f, 4), "\n")
cat("Linf estimate =", round(Linf_f, 2), "\n")
cat("Linf 95% CI   =", round(lower_Linf_f, 2), "-", round(upper_Linf_f, 2), "\n")

# ----------------------------------------
# Display numeric results for male
# ----------------------------------------
cat("\nMale estimates:\n")
cat("k estimate    =", round(estimate_k_m, 4), "\n")
cat("k 95% CI      =", round(lower_k_m, 4), "-", round(upper_k_m, 4), "\n")
cat("Linf estimate =", round(Linf_m, 2), "\n")
cat("Linf 95% CI   =", round(lower_Linf_m, 2), "-", round(upper_Linf_m, 2), "\n")
cat("========================\n")

