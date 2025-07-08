# -------------------------------------
# Function definitions for IL model
# -------------------------------------

# Inverse Gaussian density (log-scale) for increments
ff_IL <- function(Linf, k, lambda, miu, sigma, Xi, Ti, Li) {
  gLinf <- 0
  ni <- length(Xi)
  for (j in 1:ni) {
    # Calculate expected increment and IG parameters
    mu_ij <- (Linf - Li[j]) * (1 - exp(-k * Ti[j]))
    lambda_ij <- lambda
    x <- Xi[j]
    # Inverse Gaussian log-density
    llik <- 0.5*log(lambda_ij/(2*pi*x^3)) - (lambda_ij*(x - mu_ij)^2)/(2*mu_ij^2*x)
    # Lognormal prior for Linf
    lprior <- -log(sigma) - 0.5*log(2*pi) - ((log(Linf) - miu)^2)/(2*sigma^2)
    # Jacobian for transformation
    ljacob <- -log(Linf - Li[j])
    gLinf <- gLinf + llik + lprior + ljacob
  }
  return(gLinf)
}

gg_IL <- function(Linf, k, lambda, miu, sigma, Xi, Ti, Li, MinLinf, MaxLinf) {
  res <- ff_IL(Linf, k, lambda, miu, sigma, Xi, Ti, Li)
  Mi <- optimize(ff_IL, c(MinLinf, MaxLinf), k = k, lambda = lambda, miu = miu, sigma = sigma,
                 Xi = Xi, Ti = Ti, Li = Li, maximum = TRUE)$objective
  return(exp(res - Mi))
}

LL_IL <- function(theta, dat, MinLinf, MaxLinf) {
  k <- theta[1]
  lambda <- theta[2]
  miu <- theta[3]
  sigma <- theta[4]
  all.LL <- 0
  for (iid in unique(dat$LOBSTER)) {
    dati <- dat[dat$LOBSTER == iid, ]
    Xi <- dati$INC
    Li <- dati$PL
    Ti <- dati$INT / 365.25
    Mi <- optimize(ff_IL, c(MinLinf, MaxLinf), k = k, lambda = lambda, miu = miu, sigma = sigma,
                   Xi = Xi, Ti = Ti, Li = Li, maximum = TRUE)$objective
    res.int <- integrate(gg_IL, lower = MinLinf, upper = MaxLinf, k = k, lambda = lambda, miu = miu,
                         sigma = sigma, Xi = Xi, Ti = Ti, Li = Li,
                         MinLinf = MinLinf, MaxLinf = MaxLinf)$value
    res.int <- log(res.int)
    all.LL <- all.LL + res.int + Mi
  }
  return(-all.LL)
}

# -------------------------------------
# Setup shared parameters
# -------------------------------------
MinLinf <- 160
MaxLinf <- 200
init_theta_IL <- c(0.3, 4, 1, 2.88) # k, lambda, miu, sigma

# -------------------------------------
# Female
# -------------------------------------
dat <- lobster
dat_female <- dat[dat$SEX == 1, ]

res_female_IL <- optim(init_theta_IL, LL_IL, dat = dat_female,
                       MinLinf = MinLinf, MaxLinf = MaxLinf,
                       lower = c(0.01, 0.1, 0.1, 0.1),
                       upper = c(1, 100, 10, 10),
                       method = "L-BFGS-B",
                       control = list(maxit = 1000, trace = 1))

val_female_IL <- res_female_IL$par
neg_logL_female_IL <- res_female_IL$value
logL_female_IL <- -neg_logL_female_IL
num_params_IL <- 4
AIC_female_IL <- 2 * (-logL_female_IL) + 2 * num_params_IL
Linf_female_IL <- exp(val_female_IL[3] + (val_female_IL[4]^2) / 2)

cat("\n========== Female IL Model Results ==========\n")
cat("Convergence:", res_female_IL$convergence == 0, "\n")
cat("Message:", res_female_IL$message, "\n")
cat("Log-likelihood:", logL_female_IL, "\n")
cat("AIC:", AIC_female_IL, "\n")
cat("Estimated Linf:", Linf_female_IL, "\n")
cat("Parameters: k =", val_female_IL[1], ", lambda =", val_female_IL[2], 
    ", miu =", val_female_IL[3], ", sigma =", val_female_IL[4], "\n")
cat("=============================================\n")

# -------------------------------------
# Male
# -------------------------------------
dat_male <- dat[dat$SEX == 2, ]

res_male_IL <- optim(init_theta_IL, LL_IL, dat = dat_male,
                     MinLinf = MinLinf, MaxLinf = MaxLinf,
                     lower = c(0.01, 0.1, 0.1, 0.1),
                     upper = c(1, 100, 10, 10),
                     method = "L-BFGS-B",
                     control = list(maxit = 1000, trace = 1))

val_male_IL <- res_male_IL$par
neg_logL_male_IL <- res_male_IL$value
logL_male_IL <- -neg_logL_male_IL
AIC_male_IL <- 2 * (-logL_male_IL) + 2 * num_params_IL
Linf_male_IL <- exp(val_male_IL[3] + (val_male_IL[4]^2) / 2)

cat("\n========== Male IL Model Results ==========\n")
cat("Convergence:", res_male_IL$convergence == 0, "\n")
cat("Message:", res_male_IL$message, "\n")
cat("Log-likelihood:", logL_male_IL, "\n")
cat("AIC:", AIC_male_IL, "\n")
cat("Estimated Linf:", Linf_male_IL, "\n")
cat("Parameters: k =", val_male_IL[1], ", lambda =", val_male_IL[2], 
    ", miu =", val_male_IL[3], ", sigma =", val_male_IL[4], "\n")
cat("===========================================\n")
