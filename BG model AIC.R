# ---------------------------
# BG Model (Beta increments, Gamma Linf prior)
# ---------------------------

# Function for log-likelihood integrand (with Gamma prior for Linf)
ff_BG <- function(Linf, k, zeta, alpha, beta, Xi, Ti, Li) {
  gLinf <- 0
  ni <- length(Xi)
  for (j in 1:ni) {
    a <- (1 - exp(-k * Ti[j])) * (zeta - 1)
    b <- exp(-k * Ti[j]) * (zeta - 1)
    x <- Xi[j] / (Linf - Li[j])
    # Beta likelihood
    log_beta <- lgamma(a + b) + (a - 1) * log(x) + (b - 1) * log(1 - x) - lgamma(a) - lgamma(b) - log(Linf - Li[j])
    # Gamma prior for Linf
    log_prior <- (alpha - 1) * log(Linf) - Linf / beta - alpha * log(beta) - lgamma(alpha)
    gLinf <- gLinf + log_beta + log_prior
  }
  return(gLinf)
}

# Integrand for numerical integration (with normalization)
gg_BG <- function(Linf, k, zeta, alpha, beta, Xi, Ti, Li, MinLinf, MaxLinf) {
  res <- ff_BG(Linf, k, zeta, alpha, beta, Xi, Ti, Li)
  Mi <- optimize(ff_BG, c(MinLinf, MaxLinf), k = k, zeta = zeta, alpha = alpha, beta = beta,
                 Xi = Xi, Ti = Ti, Li = Li, maximum = TRUE)$objective
  return(exp(res - Mi))
}

# Negative log-likelihood for all individuals
LL_BG <- function(theta, dat, MinLinf, MaxLinf) {
  k <- theta[1]
  zeta <- theta[2]
  alpha <- theta[3]  # shape parameter for Gamma prior
  beta <- theta[4]   # scale parameter for Gamma prior
  all.LL <- 0
  for (iid in unique(dat$LOBSTER)) {
    dati <- dat[dat$LOBSTER == iid, ]
    Xi <- dati$INC
    Li <- dati$PL
    Ti <- dati$INT / 365.25
    Mi <- optimize(ff_BG, c(MinLinf, MaxLinf), k = k, zeta = zeta, alpha = alpha, beta = beta,
                   Xi = Xi, Ti = Ti, Li = Li, maximum = TRUE)$objective
    res.int <- integrate(gg_BG, lower = MinLinf, upper = MaxLinf, k = k, zeta = zeta, alpha = alpha, beta = beta,
                         Xi = Xi, Ti = Ti, Li = Li, MinLinf = MinLinf, MaxLinf = MaxLinf)$value
    res.int <- log(res.int)
    all.LL <- all.LL + res.int + Mi
  }
  return(-all.LL)
}

# ---------------------------
# Setup shared parameters
# ---------------------------
MinLinf <- 160
MaxLinf <- 200
init_theta_BG <- c(0.1, 4, 175, 10) # k, zeta, alpha, beta

# ---------------------------
# Female
# ---------------------------
dat <- lobster
dat_female <- dat[dat$SEX == 1, ]

res_female_BG <- optim(init_theta_BG, LL_BG, dat = dat_female,
                       MinLinf = MinLinf, MaxLinf = MaxLinf,
                       lower = c(0.01, 0.1, 0.1, 0.1),
                       upper = c(1, 100, 1000, 1000),
                       method = "L-BFGS-B",
                       control = list(maxit = 1000, trace = 1))

val_female_BG <- res_female_BG$par
neg_logL_female_BG <- res_female_BG$value
logL_female_BG <- -neg_logL_female_BG
num_params_BG <- 4
AIC_female_BG <- 2 * num_params_BG + 2 * neg_logL_female_BG
# For Gamma: mean = alpha * beta
Linf_female_BG <- val_female_BG[3] * val_female_BG[4]

cat("\n========== Female BG Model Results ==========\n")
cat("Convergence:", res_female_BG$convergence == 0, "\n")
cat("Message:", res_female_BG$message, "\n")
cat("Log-likelihood:", logL_female_BG, "\n")
cat("AIC:", AIC_female_BG, "\n")
cat("Estimated Linf (mean):", Linf_female_BG, "\n")
cat("Parameters: k =", val_female_BG[1], ", zeta =", val_female_BG[2], 
    ", alpha =", val_female_BG[3], ", beta =", val_female_BG[4], "\n")
cat("=============================================\n")

# ---------------------------
# Male
# ---------------------------
dat <- lobster
dat_male <- dat[dat$SEX == 2, ]
init_theta_BG_male <- c(0.1, 4, 10, 17)  # k, zeta, alpha, beta => Linf ≈ 170

res_male_BG <- optim(
  init_theta_BG_male, LL_BG, dat = dat_male,
  MinLinf = MinLinf, MaxLinf = MaxLinf,
  lower = c(0.01, 0.1, 1, 1),
  upper = c(1, 100, 500, 50),
  method = "L-BFGS-B",
  control = list(maxit = 5000, trace = 1)
)

val_male_BG <- res_male_BG$par
neg_logL_male_BG <- res_male_BG$value
logL_male_BG <- -neg_logL_male_BG
AIC_male_BG <- 2 * num_params_BG + 2 * neg_logL_male_BG
Linf_male_BG <- val_male_BG[3] * val_male_BG[4]

cat("\n========== Male BG Model Results ==========\n")
cat("Convergence:", res_male_BG$convergence == 0, "\n")
cat("Message:", res_male_BG$message, "\n")
cat("Log-likelihood:", logL_male_BG, "\n")
cat("AIC:", AIC_male_BG, "\n")
cat("Estimated Linf (mean):", Linf_male_BG, "\n")
cat("Parameters: k =", val_male_BG[1], ", zeta =", val_male_BG[2], 
    ", alpha =", val_male_BG[3], ", beta =", val_male_BG[4], "\n")
cat("===========================================\n")

