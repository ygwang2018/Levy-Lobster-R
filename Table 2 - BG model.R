# BG Model (Beta increments, Gamma Linf prior)

dat<-lobster

ff_BG <- function(Linf, k, zeta, alpha, beta, Xi, Ti, Li) {
  
  # Prior: Gamma(shape = alpha, scale = beta)
  if (Linf <= max(Li) + 1e-8) return(-1e12)
  
  log_prior <- (alpha - 1)*log(Linf) - (Linf/beta) -
    alpha*log(beta) - lgamma(alpha)
  
  out <- log_prior     # prior only ONCE
  
  for (j in seq_along(Xi)) {
    
    denom <- Linf - Li[j]
    if (denom <= 1e-12) return(-1e12)
    
    x <- Xi[j] / denom
    if (x <= 1e-12 || x >= 1 - 1e-12) return(-1e12)
    
    a <- (1 - exp(-k * Ti[j])) * (zeta - 1)
    b <- exp(-k * Ti[j]) * (zeta - 1)
    if (a <= 0 || b <= 0) return(-1e12)
    
    out <- out +
      dbeta(x, a, b, log = TRUE) -
      log(denom)
  }
  
  return(out)
}

# Normalized kernel for integration
kernel_BG <- function(Linf, Mi, k, zeta, alpha, beta, Xi, Ti, Li) {
  val <- ff_BG(Linf, k, zeta, alpha, beta, Xi, Ti, Li)
  exp(val - Mi)
}

# Integration grid
safe_integrate <- function(fun, lower, upper, n=1200) {
  xs <- seq(lower, upper, length.out=n)
  vals <- sapply(xs, fun)
  dx <- (upper - lower)/(n-1)
  sum(vals) * dx
}

# Negative log-likelihood
LL_BG <- function(theta, dat, MinLinf, MaxLinf) {
  
  k     <- theta[1]
  zeta  <- theta[2]
  alpha <- theta[3]
  beta  <- theta[4]
  
  if (k <= 0 || zeta <= 1.01 || alpha <= 0 || beta <= 0)
    return(1e12)
  
  LLtot <- 0
  
  for (id in unique(dat$LOBSTER)) {
    
    dati <- subset(dat, LOBSTER == id)
    
    Xi <- dati$INC
    Li <- dati$PL
    Ti <- dati$INT / 365.25
    
    # mode of joint density
    Mi <- optimize(
      function(Ls) ff_BG(Ls, k, zeta, alpha, beta, Xi, Ti, Li),
      interval = c(MinLinf, MaxLinf),
      maximum = TRUE
    )$objective
    
    val <- safe_integrate(
      function(Ls)
        kernel_BG(Ls, Mi, k, zeta, alpha, beta, Xi, Ti, Li),
      MinLinf, MaxLinf
    )
    
    if (!is.finite(val) || val <= 0)
      return(1e12)
    
    LLtot <- LLtot + log(val) + Mi
  }
  
  return(-LLtot)
}

dat_male <- dat[dat$SEX == 1, ]
MinLinf <- 180
MaxLinf <- 200

init_theta_BG <- c(k=0.1, zeta=4, alpha=170, beta=1.0)  # mean(L∞)=170

res_female_BG <- optim(
  par     = init_theta_BG,
  fn      = LL_BG,
  dat     = dat_female,
  MinLinf = MinLinf,
  MaxLinf = MaxLinf,
  method  = "L-BFGS-B",
  lower   = c(0.01, 1.02, 1e-4, 1e-4),
  upper   = c(2,    200,  500,  20),
  control = list(maxit = 1500, trace = 1)
)

val_female_BG <- res_female_BG$par
negLL         <- res_female_BG$value
p             <- length(val_female_BG)
AIC_female_BG <- 2*negLL + 2*p

Linf_female_BG <- val_female_BG[3] * val_female_BG[4]

print(" Female ")
print(val_female_BG)
print("AIC:", AIC_female_BG)
print("E[Linf] =", Linf_female_BG)


# Male
dat <- lobster
dat_male <- dat[dat$SEX == 2, ]
MinLinf <- 180
MaxLinf <- 200
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
p<-length(val_male_BG)
logLik <- -res_male_BG$value
AIC_male_BG <- (2 * p) - (2 * logLik)
Linf_male_BG <- val_male_BG[3] * val_male_BG[4]
AIC_male_BG
Linf_male_BG
val_male_BG[1]

#=== Male  ======
print("Convergence:", res_male_BG$convergence == 0)
print("Log-likelihood:", logL_male_BG)
print("AIC:", AIC_male_BG)
print("Estimated Linf (mean):", Linf_male_BG)
print("Parameters: k =", val_male_BG[1], ", zeta =", val_male_BG[2], 
    ", alpha =", val_male_BG[3], ", beta =", val_male_BG[4])

