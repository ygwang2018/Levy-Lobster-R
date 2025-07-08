# ---------------------------
# Function Definitions
# ---------------------------

ff_GL <- function(Linf, k, alpha, miu, sigma, Xi, Ti, Li) {
  gLinf <- 0
  ni <- length(Xi)
  for (j in 1:ni) {
    # Gamma increment likelihood
    mu <- (Linf - Li[j]) * (1 - exp(-k * Ti[j]))
    shape <- alpha
    scale <- mu / alpha
    log_inc <- dgamma(Xi[j], shape=shape, scale=scale, log=TRUE)
    # Lognormal prior for Linf
    log_prior <- dlnorm(Linf, meanlog=miu, sdlog=sigma, log=TRUE)
    gLinf <- gLinf + log_inc + log_prior
  }
  return(gLinf)
}

gg_GL <- function(Linf, k, alpha, miu, sigma, Xi, Ti, Li, MinLinf, MaxLinf) {
  res <- ff_GL(Linf, k, alpha, miu, sigma, Xi, Ti, Li)
  Mi <- optimize(ff_GL, c(MinLinf, MaxLinf), k=k, alpha=alpha, miu=miu, sigma=sigma, Xi=Xi, Ti=Ti, Li=Li, maximum=TRUE)$objective
  return(exp(res - Mi))
}

LL_GL <- function(theta, dat, MinLinf, MaxLinf) {
  k <- theta[1]
  alpha <- theta[2]
  miu <- theta[3]
  sigma <- theta[4]
  all.LL <- 0
  for (iid in unique(dat$LOBSTER)) {
    dati <- dat[dat$LOBSTER == iid, ]
    Xi <- dati$INC
    Li <- dati$PL
    Ti <- dati$INT / 365.25
    # Skip increments that cannot be explained by any Linf
    if (any(Li + Xi >= MaxLinf)) next
    Mi <- optimize(ff_GL, c(MinLinf, MaxLinf), k=k, alpha=alpha, miu=miu, sigma=sigma, Xi=Xi, Ti=Ti, Li=Li, maximum=TRUE)$objective
    res.int <- tryCatch({
      integrate(gg_GL, lower=MinLinf, upper=MaxLinf, k=k, alpha=alpha, miu=miu, sigma=sigma,
                Xi=Xi, Ti=Ti, Li=Li, MinLinf=MinLinf, MaxLinf=MaxLinf, rel.tol=1e-4)$value
    }, error=function(e) NA)
    if (is.na(res.int) || res.int <= 0) next
    res.int <- log(res.int)
    all.LL <- all.LL + res.int + Mi
  }
  return(-all.LL)
}

# ---------------------------
# Model Setup
# ---------------------------
MinLinf <- 160
MaxLinf <- 200
init_theta <- c(0.1, 4, 5, 0.5) # k, alpha, miu, sigma (set sensible starting values)
num_params <- 4

# ---------------------------
# Female
# ---------------------------
dat_female <- lobster[lobster$SEX == 1, ]
res_female_GL <- optim(init_theta, LL_GL, dat=dat_female,
                       MinLinf=MinLinf, MaxLinf=MaxLinf,
                       lower=c(0.01, 0.1, 0.1, 0.01),
                       upper=c(1, 100, 10, 10),
                       method="L-BFGS-B",
                       control=list(maxit=1000, trace=1))
val_female_GL <- res_female_GL$par
neg_logL_female_GL <- res_female_GL$value
logL_female_GL <- -neg_logL_female_GL
AIC_female_GL <- 2 * num_params + 2 * neg_logL_female_GL
Linf_female_GL <- exp(val_female_GL[3] + (val_female_GL[4]^2) / 2)

cat("\n========== Female GL Results ==========\n")
cat("Convergence:", res_female_GL$convergence == 0, "\n")
cat("Message:", res_female_GL$message, "\n")
cat("Log-likelihood:", logL_female_GL, "\n")
cat("AIC:", AIC_female_GL, "\n")
cat("Estimated Linf:", Linf_female_GL, "\n")
cat("Parameters: k =", val_female_GL[1], ", alpha =", val_female_GL[2], 
    ", miu =", val_female_GL[3], ", sigma =", val_female_GL[4], "\n")
cat("=======================================\n")

# ---------------------------
# Male
# ---------------------------
dat_male <- lobster[lobster$SEX == 2, ]
res_male_GL <- optim(init_theta, LL_GL, dat=dat_male,
                     MinLinf=MinLinf, MaxLinf=MaxLinf,
                     lower=c(0.01, 0.1, 0.1, 0.01),
                     upper=c(1, 100, 10, 10),
                     method="L-BFGS-B",
                     control=list(maxit=1000, trace=1))
val_male_GL <- res_male_GL$par
neg_logL_male_GL <- res_male_GL$value
logL_male_GL <- -neg_logL_male_GL
AIC_male_GL <- 2 * num_params + 2 * neg_logL_male_GL
Linf_male_GL <- exp(val_male_GL[3] + (val_male_GL[4]^2) / 2)

cat("\n========== Male GL Results ==========\n")
cat("Convergence:", res_male_GL$convergence == 0, "\n")
cat("Message:", res_male_GL$message, "\n")
cat("Log-likelihood:", logL_male_GL, "\n")
cat("AIC:", AIC_male_GL, "\n")
cat("Estimated Linf:", Linf_male_GL, "\n")
cat("Parameters: k =", val_male_GL[1], ", alpha =", val_male_GL[2], 
    ", miu =", val_male_GL[3], ", sigma =", val_male_GL[4], "\n")
cat("=====================================\n")
