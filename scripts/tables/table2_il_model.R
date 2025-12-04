# Log-integrand for one individual
dat<-lobster
ff_IL <- function(Linf, k, lambda, mu_L, sigma_L, Xi, Ti, Li) {
  
  if (Linf <= max(Li) + 1e-10) return(-1e12)
  
  logprior <- -log(sigma_L) - log(Linf) - 0.5*log(2*pi) -
    ((log(Linf) - mu_L)^2) / (2*sigma_L^2)
  
  ll <- logprior
  
  for (j in seq_along(Xi)) {
    
    # IG mean
    mu_ij <- (Linf - Li[j]) * (1 - exp(-k * Ti[j]))
    if (mu_ij <= 0) return(-1e12)
    
    x <- Xi[j]
    if (x <= 0) return(-1e12)
    
    # Inverse Gaussian log-density
    log_igg <- 0.5 * log(lambda/(2*pi*x^3)) -
      lambda * (x - mu_ij)^2 / (2 * mu_ij^2 * x)
    
    # Jacobian from Δ = Λ (L∞ - L)
    ljacob <- log(Linf - Li[j])
    
    ll <- ll + log_igg + ljacob
  }
  
  return(ll)
}

# Normalised kernel
kernel_IL <- function(Linf, Mi, k, lambda, mu_L, sigma_L, Xi, Ti, Li) {
  val <- ff_IL(Linf, k, lambda, mu_L, sigma_L, Xi, Ti, Li)
  exp(val - Mi)
}

# Grid integration
safe_integrate <- function(fun, lower, upper, n = 1200) {
  xs <- seq(lower, upper, length.out=n)
  vals <- sapply(xs, fun)
  dx <- (upper - lower)/(n - 1)
  sum(vals) * dx
}

LL_IL <- function(theta, dat, MinLinf, MaxLinf) {
  
  k       <- theta[1]
  lambda  <- theta[2]
  mu_L    <- theta[3]
  sigma_L <- theta[4]
  
  if (k <= 0 || lambda <= 0 || sigma_L <= 0)
    return(1e12)
  
  LLtot <- 0
  
  for (id in unique(dat$LOBSTER)) {
    
    dati <- subset(dat, LOBSTER == id)
    Xi <- dati$INC
    Li <- dati$PL
    Ti <- dati$INT / 365.25
    
    Mi <- optimize(
      function(Ls) ff_IL(Ls, k, lambda, mu_L, sigma_L, Xi, Ti, Li),
      interval = c(MinLinf, MaxLinf),
      maximum = TRUE
    )$objective
    
    val <- safe_integrate(
      function(Ls)
        kernel_IL(Ls, Mi, k, lambda, mu_L, sigma_L, Xi, Ti, Li),
      MinLinf, MaxLinf
    )
    
    if (!is.finite(val) || val <= 0)
      return(1e12)
    
    LLtot <- LLtot + log(val) + Mi
  }
  
  return(-LLtot)
}

# FEMALE FIT

dat_female <- subset(lobster, SEX==1)

MinLinf <- 100
MaxLinf <- 200

init_theta_IL <- c(
  k = 0.1,
  lambda = 2,
  mu_L = log(180),
  sigma_L = 0.15
)

res_female_IL <- optim(
  par     = init_theta_IL,
  fn      = LL_IL,
  dat     = dat_female,
  MinLinf = MinLinf,
  MaxLinf = MaxLinf,
  method  = "L-BFGS-B",
  lower   = c(0.1, 0.1, log(120), 0.05),
  upper   = c(0.5, 50,  log(260), 0.50),
  control = list(maxit = 3000, trace = 1)
)

val_female_IL <- res_female_IL$par
negLL_female  <- res_female_IL$value
p_female      <- length(val_female_IL)

AIC_female_IL <- -2*negLL_female + 2*p_female
Linf_female_IL <- exp(val_female_IL[3] + val_female_IL[4]^2/2)

#  FEMALE  
print(val_female_IL)
print("AIC:", AIC_female_IL)
print("E[Linf]:", Linf_female_IL)

############ Male
dat<-lobster
dat_male <- dat[dat$SEX == 2, ]
MinLinf <- 100
MaxLinf <- 200

init_theta_IL <- c(0.1, 0.1, 0.1, 0.01)   # k, lambda, mu_L, sigma_L

res_male_IL <- optim(init_theta_IL, LL_IL, dat = dat_male,
                     MinLinf = MinLinf, MaxLinf = MaxLinf,
                     lower = c(0.1, 0.1, 0.1, 0.1),
                     upper = c(0.5, 100, 10, 10),
                     method = "L-BFGS-B",
                     control = list(maxit = 1000, trace = 1))

val_male_IL <- res_male_IL$par
neg_logL_male_IL <- res_male_IL$value
AIC_male_IL <- 2 * (-neg_logL_male_IL) + 2 * length(val_male_IL)
Linf_male_IL <- exp(val_male_IL[3] + (val_male_IL[4]^2) / 2)
AIC_male_IL
Linf_male_mean
val_male_IL[1]
val_male_IL

print(" Male Results ")
print("Convergence:", res_male_IL$convergence == 0)
print("AIC:", AIC_male_IL)
print("Estimated Linf:", Linf_male_IL)
print("Parameters: k =", val_male_IL[1], ", lambda =", val_male_IL[2], 
    ", miu =", val_male_IL[3], ", sigma =", val_male_IL[4])
